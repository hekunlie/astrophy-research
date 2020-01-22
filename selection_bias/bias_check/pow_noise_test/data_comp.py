import h5py
from mpi4py import MPI
import numpy
from sys import argv
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


src_path = "/mnt/perc/hklee/bias_check/data_from_pi/pow_noise_test"

dst_path = "%s/%s/data_%d_epsf.hdf5"%(src_path, argv[1], rank)

argc = len(argv)

data = numpy.zeros((20000000,4))

for i in range(2,argc):
    h5f = h5py.File(src_path + "/data_4/data_%d_%s_epsf.hdf5"%(rank, argv[i]))
    temp = h5f["/data"][()]
    h5f.close()

    data += temp
    if rank == 0:
        print(argv[i])

h5f = h5py.File(dst_path,"w")
h5f["/data"] = data
h5f.close()





