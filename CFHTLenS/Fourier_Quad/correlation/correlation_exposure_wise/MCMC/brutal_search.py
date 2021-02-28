import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import correlation_function_tool as cf_tool
import tool_box
import numpy
import time
import h5py
from mpi5py import MPI


t1 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

expo = int(argv[1])
expo_type = ["diff_expo","same_expo"][expo]

h = 67.4
omega_bm_num = 10
omega_bm_bin = numpy.linspace(0.01,0.1,omega_bm_num)
omega_cm_num = 40
omega_cm_bin = numpy.linspace(0.01,0.4,omega_cm_num)
sigma8_num = 71
sigma8_bin = numpy.linspace(0.3,1,sigma8_num)

chisq = numpy.zeros((omega_bm_num, omega_cm_num, sigma8_num))

ijk = []
for i in range(omega_bm_num):
    for j in range(omega_cm_num):
        for k in range(sigma8_num):
            ijk.append([i,j,k])

ijk_sub = tool_box.alloc(ijk, numprocs)[rank]


h5f = h5py.File("./data/zhist.hdf5","r")
zehist = h5f["/zhist"][()]
zebin = h5f["/zbin"][()]
zebin_cent = h5f["/zbin_cent"][()]
h5f.close()

redshift_bin_num = zehist.shape[0]
tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)


resample_num = 200

h5f = h5py.File("./data/result_cache_%d_%s.hdf5"%(resample_num,expo_type),"r")

theta = h5f["/theta"][()]
# only \xi_+
xi = h5f["/xip"][()]
cov_inv = h5f["/cov_inv"][()]
h5f.close()

# arcmin to radian & degree
theta_radian = theta/60/180*numpy.pi
theta_degree = theta/60

data_num = theta_radian.shape[0]
theta_radian = theta_radian.reshape((tomo_panel_num, int(data_num*1.0/tomo_panel_num)))
theta_degree = theta_degree.reshape((tomo_panel_num, int(data_num*1.0/tomo_panel_num)))


ell = tool_box.set_bin_log(10, 20000, 10000)

for tag in ijk_sub:
    i, j, k = tag
    omega_bm0, omega_cm0, sigma8 = omega_bm_bin[i], omega_cm_bin[j], sigma8_bin[k]

    xi_predict = cf_tool.get_tomo_xi_ccl(sigma8, omega_cm0, omega_bm0, h,
                                         theta_degree, xi, cov_inv, zebin_cent, zehist, ell)[0]

    diff = xi - xi_predict.flatten()
    chisq[i,j,k] = - 0.5 * numpy.dot(diff, numpy.dot(cov_inv, diff))

comm.Barrier()

if rank > 0:
    comm.Send([chisq, MPI.DOUBLE], dest=0, tag=rank)
else:
    for ir in range(1, numprocs):
        recv_buf = numpy.empty((omega_bm_num, omega_cm_num, sigma8_num), dtype=numpy.float64)
        comm.Recv(recv_buf, source=ir, tag=ir)

        chisq = chisq + recv_buf
    h5f = h5py.File("./data/chisq.hdf5","w")
    h5f["/chisq"] = chisq
    h5f["/omega_bm_bin"] = omega_bm_bin
    h5f["/omega_cm_bin"] = omega_cm_bin
    h5f["/sigma8"] = sigma8_bin
    h5f.close()

comm.Barrier()