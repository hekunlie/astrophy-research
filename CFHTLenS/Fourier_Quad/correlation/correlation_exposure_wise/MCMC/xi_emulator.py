import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import correlation_function_tool as cf_tool
import tool_box
import numpy
import time
import h5py
from mpi4py import MPI



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

st = time.time()

redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)

h5f = h5py.File("./data/zhist.hdf5","r")
zehist = h5f["/zhist"][()]
zebin = h5f["/zbin"][()]
zebin_cent = h5f["/zbin_cent"][()]
h5f.close()

redshift_bin_num = 6
tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)

h = 0.674

ell = tool_box.set_bin_log(10, 20000, 10000)

used_zbins = numpy.array([0,1,1,1,1,1])

sigma8_num = 250
omega_c_num = 100
omega_b_num = 50
theta_num = 15

sigma8_list = numpy.linspace(0.3, 2.5, sigma8_num)
omega_c_list = numpy.linspace(0.01, 0.8,omega_c_num)
omega_b_list = numpy.linspace(0.01,0.2,omega_b_num)
theta_deg = numpy.zeros((15, theta_num))

for i in range(15):
    theta_deg[i] = tool_box.set_bin_log(0.8, 60, theta_num)/60
theta_deg = theta_deg.flatten()

task_list = []
for i in range(sigma8_num):
    for j in range(omega_c_num):
        for k in range(omega_b_num):
            task_list.append([i,j,k])

my_task = tool_box.alloc(task_list, numprocs, "mean")[rank]

my_xip = numpy.zeros((sigma8_num, omega_c_num, omega_b_num, int(theta_num*15)), dtype=numpy.float32)
logger = tool_box.get_logger("./logs/log_%d.dat"%rank)

for tag, ijk in enumerate(my_task):

    t1 = time.time()

    i,j,k = ijk
    sigma8 = sigma8_list[i]
    omega_c = omega_c_list[j]
    omega_b = omega_b_list[k]
    logger.info("%d. %.4f %.4f %.4f." % (tag, sigma8, omega_c, omega_b))
    try:
        ccl_xip = cf_tool.get_tomo_xi_ccl_2(sigma8, omega_c, omega_b, h, zebin_cent, zehist, theta_deg, theta_num, ell, used_zbins)
        my_xip[i, j, k] = ccl_xip
        t2 = time.time()
        logger.info("%d. %.2f sec." % (tag, t2 - t1))
    except:
        my_xip[i, j, k] = - 99
        t2 = time.time()
        logger.info("%d. Failed. %.2f sec." % (tag, t2 - t1))



numpy.savez("./data/xip_%d.npz"%rank, my_xip)

ed = time.time()
comm.Barrier()

print("%d. %.2f sec"%(rank, ed-st))

