import numpy as numpy
from mpi4py import MPI
import h5py
import time

# comm it the COMMON WORLD where all the CPUs in and communicate
comm = MPI.COMM_WORLD
# the identity of the CPU, different CPU has different rank, from 0 ~ 'cpus'
rank = comm.Get_rank()
# how many CPUs you used
cpus = comm.Get_size()



#################### Send/Recv #################
# the data array in each CPU
cata_data = numpy.zeros((4,4))
# the array shape, tuple
data_sp = cata_data.shape
# gather all the shapes from each CPU
# the 'data_sps' is a list of tuple collected from each CPU
data_sps = comm.gather(data_sp, root=0)

# rank 0 will collect all data, the others send data to it.
if rank > 0:
    # !!!! remember the data type, MPI.DOUBLE, MPI.FLOAT, ...
    # or it will raise an error, Keyerror
    comm.Send([cata_data, MPI.DOUBLE], dest=0, tag=rank)
else:
    # ############ this part may not be needed #######
    if data_sp[0] > 1 and data_sp[1] > 1:
        stack_pool = [cata_data]
    else:
        stack_pool = []
    # ################################################

    # receive the data from other CPUs
    # !!!! the start points is 1 in range() not 0
    for procs in range(1, cpus):
        # prepare a buffer for the data, the shape must be the same
        # with that of what the other CPUs send, you have collected them in 'data_sps'
        recvs = numpy.empty(data_sps[procs], dtype=numpy.double)
        # receive it using the buffer,
        comm.Recv(recvs, source=procs, tag=procs)
        # then do whatever you want ...
        if data_sps[procs][0] > 1 and data_sps[procs][1] > 1:
            stack_pool.append(recvs)



#################### shared memory #################
# show how to create two contiguous shared blocks in memory
# length of double
itemsize = MPI.DOUBLE.Get_size()
element_num = 10
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
win2 = MPI.Win.Allocate_shared(nbytes*5, itemsize, comm=comm)
# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf1, itemsize = win1.Shared_query(0)
buf2, itemsize = win2.Shared_query(0)

# create a numpy array from buf
# code can run successfully without the following step
# buf = np.array(buf, dtype='float64', copy=False) # may be redundant
# "d" means double = 'float64'
ary1 = numpy.ndarray(buffer=buf1, dtype='d', shape=(element_num,)) # array filled with zero
ary2 = numpy.ndarray(buffer=buf2, dtype='d', shape=(element_num, 5))
for i in range(5):
    ary2[rank,i] = rank
# the rank 1 changes the array
if rank == 1:
    # show how to read data from a hdf5 file and
    # assign to the created block in memory

    # f = h5py.File("./para.hdf5","r")

    # ary2 = f["/data"].value does not work,
    # it must be ary2[:10] = f["/data"].value
    # or ary2[:10,0] = f["/data"].value[:10,0]
    # the shapes must match each other

    # ary2[:10] = f["/data"].value
    # f.close()

    # change the array
    ary1[:element_num] = numpy.arange(element_num, dtype='i')
# wait in process rank 0 of comm until process 1 has written to the array
comm.Barrier() # necessary

# ary_sum1 = numpy.sum(ary1**2)
# ary_sum12 = numpy.sum(ary1 ** 2)
# ary_cos = numpy.cos(ary1**2)*numpy.sin(ary1**2)
# ary_sum2 = numpy.sum(ary2)
time.sleep(0.001)
# the you will see that the changed array will be seen by each process

print(rank, ary1.shape, ary2)#, ary_sum1, ary_sum12, ary_sum2)



