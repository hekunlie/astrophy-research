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

data_path, theta_tag = argv[1], int(argv[2])

theta = numpy.linspace(0, numpy.pi, 17)[theta_tag]

rot_2theta = [[numpy.cos(2*theta), -numpy.sin(2*theta)],
              [numpy.sin(2*theta), numpy.cos(2*theta)]]

rot_4theta = [[numpy.cos(4*theta), -numpy.sin(4*theta)],
              [numpy.sin(4*theta), numpy.cos(4*theta)]]

data_nm_1 = data_path + "/data_noise_free_%d.hdf5"%rank
data_nm_2 = data_path + "/data_gal_noise_cross_term_%d.hdf5"%rank
data_nm_3 = data_path + "/data_pure_gal_noise_cross_term_est_r_%d.hdf5"%rank
# data_nm_2 = data_path + "/data_noise_free_%d.hdf5"%rank

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

h5f = h5py.File(data_nm_3, "r")
mg1_ct = h5f["/mg1"][()]
mg2_ct = h5f["/mg2"][()]
mn_ct = h5f["/mn"][()]
mu_ct = h5f["/mu"][()]
mv_ct = h5f["/mv"][()]
h5f.close()

mg1 = mg1 + mg1_ct*rot_2theta[0][0] + mg2_ct*rot_2theta[0][1]
mg2 = mg2 + mg1_ct*rot_2theta[1][0] + mg2_ct*rot_2theta[1][1]

mn = mn + mn_ct
mu = mu + mu_ct*rot_4theta[0][0] + mv_ct*rot_4theta[0][1]
mv = mv + mu_ct*rot_4theta[1][0] + mv_ct*rot_4theta[1][1]

if rank == 0:
    print(theta)
    if not os.path.exists(data_path + "/mix2"):
        os.makedirs(data_path + "/mix2")

comm.Barrier()

h5f = h5py.File(data_path + "/mix2/data_mix_%d.hdf5"%rank,"w")
h5f["/mg1"] = numpy.float32(mg1)
h5f["/mg2"] = numpy.float32(mg2)
h5f["/mn"] = numpy.float32(mn)
h5f["/mu"] = numpy.float32(mu)
h5f["/mv"] = numpy.float32(mv)
h5f.close()