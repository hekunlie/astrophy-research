import numpy as np
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()

# length of double
itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elememts
    nbytes = 10*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf, itemsize = win.Shared_query(0)


# create a numpy array from buf
# code can run successfully without the following step
# buf = np.array(buf, dtype='float64', copy=False) # may be redundant
# "d" means double = 'float64'
ary = np.ndarray(buffer=buf, dtype='d', shape=(10,)) # array filled with zero


# the rank 1 changes the array
if rank == 1:
    ary[:5] = np.arange(5, dtype='i')
# wait in process rank 0 of comm until process 1 has written to the array
comm.Barrier() # necessary

# the you will see that the changed array will be seen by each process
print(rank, ary)



