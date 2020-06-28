import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append('%s/work/mylib/' % my_home)
import h5py
from mpi4py import MPI
import tool_box
import numpy

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

data_path = argv[1]

data_nm_1 = data_path + "/data_noise_free_%d.hdf5"%rank
data_nm_2 = data_path + "/data_gal_noise_cross_term_%d.hdf5"%rank
# data_nm_2 = data_path + "/data_noise_free_%d.hdf5"%rank

set_name_1 = [["/mg1", "/mg2", "/mn", "/mu"]]

h5f = h5py.File(data_nm_1, "r")
mg1 = h5f["/mg1"][()]
mg2 = h5f["/mg2"][()]
mn = h5f["/mn"][()]
mu = h5f["/mu"][()]
mv = h5f["/mv"][()]
h5f.close()

h5f = h5py.File(data_nm_2,"r")
mg1 += h5f["/mg1"][()]
mg2 += h5f["/mg2"][()]
mn += h5f["/mn"][()]
mu += h5f["/mu"][()]
mv += h5f["/mv"][()]
h5f.close()

if rank == 0:
    if not os.path.exists(data_path + "/mix"):
        os.makedirs(data_path + "/mix")

comm.Barrier()

h5f = h5py.File(data_path + "/mix/data_mix_%d.hdf5"%rank,"w")
h5f["/mg1"] = numpy.float32(mg1)
h5f["/mg2"] = numpy.float32(mg2)
h5f["/mn"] = numpy.float32(mn)
h5f["/mu"] = numpy.float32(mu)
h5f["/mv"] = numpy.float32(mv)
h5f.close()