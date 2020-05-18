import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
import numpy
from sys import path, argv
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
path.append('%s/work/mylib' % my_home)
from Fourier_Quad import Fourier_Quad
# import h5py
# from plot_tool import Image_Plot
import tool_box
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

source_num = int(argv[1])*10000
sigma_1 = float(argv[2])
sigma_2 = float(argv[3])
signal_num = numprocs
signals = numpy.linspace(-0.05, 0.05, signal_num)

itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = 2*signal_num*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
result = numpy.ndarray(buffer=buf1, dtype='d', shape=(2, signal_num)) # array filled with zero

fq = Fourier_Quad(12,123)
n = numpy.ones((source_num, ))
# for i in range(signal_num):
source = numpy.random.normal(signals[rank], sigma_1, source_num) + numpy.random.normal(-signals[rank]/100, sigma_2, source_num)
signal_est = fq.find_shear(source, n, 8,scale=100, left=-0.08, right=0.08)[:2]
result[:, rank] = signal_est
print(rank, signal_est)
comm.Barrier()
if rank == 0:
    # result[2] = signals
    print(signals)
    print(result)
    mc = numpy.array(tool_box.data_fit(signals, result[0], result[1]))
    mc[0] = mc[0] - 1
    print(mc)
# img = Image_Plot()
# img.subplots(1,1)
# img.axs[0][0].errorbar(signals, result[0], result[1])
# img.axs[0][0].plot([-0.06,0.06],[-0.06, 0.06])
# img.show_img()