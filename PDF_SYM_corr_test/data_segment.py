import numpy
import h5py
from sys import argv
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

if rank == 0:
    data_type = "noise_free"
else:
    data_type = "noisy_cpp"
total_path = "/mnt/perc/hklee/PDF_test"
for i in range(30):
    h5f = h5py.File(total_path + "/data_ori/data_%d_%s.hdf5"%(i,data_type),"r")
    src = h5f["/data"][()][:,:4]

    h5f.close()
    h5f = h5py.File(total_path + "/data/data_%d_%s.hdf5"%(i,data_type),"w")
    for j in range(4):
        st,ed = int(j*25000000), int((j+1)*25000000)
        h5f["/data_%d"%j] = src[st:ed]
    h5f.close()
    print(i,src.shape)