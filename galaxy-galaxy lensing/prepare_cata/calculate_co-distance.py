import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
path.append('%s/work/mylib/' % my_home)
from astropy.cosmology import FlatLambdaCDM
import numpy
import h5py
import time
from mpi4py import MPI
import tool_box
from plot_tool import Image_Plot

""" calculate comoving distance"""

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


# data path
data_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/redshift.hdf5"

# parameters
omega_m0 = 0.31
omega_lambda0 = 1 - omega_m0
h = 0.7
C_0_hat = 2.99792458
H_0 = 100*h
coeff = 1000*C_0_hat

# total point number
total_z = 2000001
dz = 5
# redshift array
redshift = numpy.zeros((total_z,))
for i in range(total_z):
    redshift[i] = i*dz/1000000.

# task distribution
total_task = numpy.arange(0,total_z)
my_task = tool_box.allot(total_task, cpus)[rank]


# shared buffer
itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    nbytes = total_z*itemsize*2
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
# comoving distance (Mpc/h), integrate part of co-distance
data_buf = numpy.ndarray(buffer=buf1, dtype='d', shape=(total_z, 2))
comm.Barrier()

# calculate
cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)
t1 = time.time()
for i in my_task:
    # Mpc
    com_dist = cosmos.comoving_distance(redshift[i])
    # Mpc/h
    data_buf[i,0] = com_dist.value*h
    # the integrate in the distance
    data_buf[i,1] = com_dist.value*h/coeff

    if rank == 0:
        print(redshift[i], com_dist, data_buf[i,0], data_buf[i,1])
t2 = time.time()
print("%d, %.2f sec"%(rank,t2-t1))

comm.Barrier()
if rank == 0:
    h5f = h5py.File(data_path, "w")
    h5f["/OM0_H0_C"] = numpy.array([omega_m0, H_0, C_0_hat])
    h5f["/Z"] = redshift
    h5f["/DISTANCE"] = data_buf[:,0]
    h5f["/DISTANCE_INTEG"] = data_buf[:,1]
    h5f.close()

    img = Image_Plot()
    img.subplots(1,2)
    img.axs[0][0].plot(redshift,data_buf[:,0])
    img.axs[0][1].plot(redshift,data_buf[:,1])
    img.save_img("test.png")
    img.close_img()
