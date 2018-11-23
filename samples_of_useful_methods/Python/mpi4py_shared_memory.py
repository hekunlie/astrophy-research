import numpy as np
from mpi4py import MPI
import h5py

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# show how to create two contiguous shared blocks in memory
# length of double
itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = 10*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
win2 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf1, itemsize = win1.Shared_query(0)
buf2, itemsize = win2.Shared_query(0)

# create a numpy array from buf
# code can run successfully without the following step
# buf = np.array(buf, dtype='float64', copy=False) # may be redundant
# "d" means double = 'float64'
ary1 = np.ndarray(buffer=buf1, dtype='d', shape=(10,)) # array filled with zero
ary2 = np.ndarray(buffer=buf2, dtype='d', shape=(10,1))

# the rank 1 changes the array
if rank == 1:
    # show how to read data from a hdf5 file and
    # assign to the created block in memory
    f = h5py.File("./para.hdf5","r")
    # ary2 = f["/data"].value does not work,
    # it must be ary2[:10] = f["/data"].value
    # or ary2[:10,0] = f["/data"].value[:10,0]
    # the shapes must match each other
    ary2[:10] = f["/data"].value
    f.close()

    # change the array
    ary1[:5] = np.arange(5, dtype='i')

# wait in process rank 0 of comm until process 1 has written to the array
comm.Barrier() # necessary

# the you will see that the changed array will be seen by each process
print(rank, ary1, ary2)



