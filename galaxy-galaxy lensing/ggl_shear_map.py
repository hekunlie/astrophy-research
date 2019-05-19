import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import time
from Fourier_Quad import Fourier_Quad
import tool_box
import plot_tool
import h5py
from mpi4py import MPI
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

area_id = int(argv[1])
foreground = argv[2]

itemsize = MPI.DOUBLE.Get_size()
element_num = 13*4
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf1, itemsize = win1.Shared_query(0)


# create a numpy array from buf
# code can run successfully without the following step
# buf = np.array(buf, dtype='float64', copy=False) # may be redundant
# "d" means double = 'float64'
result = numpy.ndarray(buffer=buf1, dtype='d', shape=(4,13)) # array filled with zero