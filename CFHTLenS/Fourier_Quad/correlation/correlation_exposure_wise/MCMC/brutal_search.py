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
import warnings

warnings.filterwarnings('error')

t1 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

data_path = argv[1]
expo = int(argv[2])
expo_type = ["diff_expo","same_expo"][expo]

h = 0.674
omega_bm_num = 10
omega_bm_bin = numpy.linspace(0.01,0.1,omega_bm_num)
omega_cm_num = 91
omega_cm_bin = numpy.linspace(0.1,1,omega_cm_num)
sigma8_num = 71
sigma8_bin = numpy.linspace(0.3, 1,sigma8_num)

chisqs = numpy.zeros((omega_bm_num, omega_cm_num, sigma8_num))

ijk = []
for i in range(omega_bm_num):
    for j in range(omega_cm_num):
        for k in range(sigma8_num):
            ijk.append([i,j,k])

ijk_sub = tool_box.alloc(ijk, numprocs)[rank]


h5f = h5py.File("./%s/zhist.hdf5"%data_path,"r")
zehist = h5f["/zhist"][()]
zebin = h5f["/zbin"][()]
zebin_cent = h5f["/zbin_cent"][()]
h5f.close()

redshift_bin_num = zehist.shape[0]
tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)


resample_num = 200

h5f = h5py.File("./%s/result_cache_%d_%s.hdf5"%(data_path,resample_num, expo_type), "r")
theta = h5f["/theta"][()]
# only \xi_+
xi = h5f["/xi_p"][()]
cov_inv = h5f["/inv_cov_xi_p"][()]
used_zbins = h5f["/used_zbin"][()]
h5f.close()

print("Used Z bins: ", used_zbins)


# arcmin to radian & degree
theta_radian = theta/60/180*numpy.pi
theta_degree = theta/60

theta_num_per_bin = 5

data_num = theta_radian.shape[0]

print("Data vec: ", theta_radian.shape)
# print(theta)
# print(theta_radian)
# print(theta_degree)

logger = tool_box.get_logger("./logs/log_%d.dat"%rank)

ell = tool_box.set_bin_log(5, 15000, 200)

# xi_cache = numpy.zeros((len(ijk_sub), data_num)) - 1
count = 0
total_step = len(ijk_sub)
for tag in ijk_sub:
    i, j, k = tag
    omega_bm0, omega_cm0, sigma8 = omega_bm_bin[i], omega_cm_bin[j], sigma8_bin[k]
    try:
        xi_predict = cf_tool.get_tomo_xi_ccl_2(sigma8, omega_cm0, omega_bm0, h, zebin_cent, zehist,
                                               theta_degree, theta_num_per_bin, ell, used_zbins)

        diff = xi - xi_predict.flatten()
        chi_sq = - 0.5 * numpy.dot(diff, numpy.dot(cov_inv, diff))

        # xi_cache[count] = xi_predict
        logger.info("%dth/%d step. \Xi^2 %.4f Omega_b =%.3f,Omega_c =%.3f, sigma8 =%.3f" %
                    (count+1, total_step, chi_sq, omega_bm0, omega_cm0, sigma8))
    except:
        logger.info("%dth/%d step. Failed at Omega_b =%.3f,Omega_c =%.3f, sigma8 =%.3f"%(count+1,total_step, omega_bm0, omega_cm0, sigma8))
        print(omega_bm0, omega_cm0, sigma8)
        chi_sq = - numpy.inf

    chisqs[i, j, k] = chi_sq
    count += 1

# numpy.savez("./xi_pts/xi_%d.npz"%rank, xi_cache, ijk_sub)

comm.Barrier()

if rank > 0:
    comm.Send([chisqs, MPI.DOUBLE], dest=0, tag=rank)
else:
    for ir in range(1, numprocs):
        recv_buf = numpy.empty((omega_bm_num, omega_cm_num, sigma8_num), dtype=numpy.float64)
        comm.Recv(recv_buf, source=ir, tag=ir)

        chisqs = chisqs + recv_buf

    h5f = h5py.File("./%s/chisq_%s.hdf5"%(data_path,expo_type),"w")
    h5f["/chisq"] = chisqs
    h5f["/omega_bm_bin"] = omega_bm_bin
    h5f["/omega_cm_bin"] = omega_cm_bin
    h5f["/sigma8"] = sigma8_bin
    h5f.close()

comm.Barrier()