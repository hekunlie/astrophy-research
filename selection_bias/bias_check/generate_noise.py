import numpy
import h5py
from mpi4py import MPI
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

parent_path = "/mnt/ddnfs/data_users/hkli/noise_img/%d"%rank
noise_sigma = 60
stamp_size = 64
stamp_nx,stamp_ny = 100,100
stamp_num = stamp_nx*stamp_ny
noise_pts_on_chip = stamp_num*stamp_size*stamp_size
shear_num = 20
total_chip_num = 1000

for i in range(total_chip_num):
    seed1 = (rank+1)*20000 + i*3
    seed2 = seed1 + 10000
    rng1 = numpy.random.RandomState(seed1)
    rng2 = numpy.random.RandomState(seed2)

    noise_1 = rng1.normal(0,noise_sigma,noise_pts_on_chip).reshape((stamp_size*stamp_ny, stamp_size*stamp_ny))
    noise_2 = rng2.normal(0,noise_sigma,noise_pts_on_chip).reshape((stamp_size*stamp_ny, stamp_size*stamp_ny))

    h5f = h5py.File(parent_path + "/chip_%d.hdf5"%i,"w")
    h5f["/noise_1"] = noise_1
    h5f["/noise_2"] = noise_2
    h5f.close()