import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import tool_box
import correlation_function_tool as cf_tool
import time
import numpy
import emcee
import time
import h5py
from mpi4py import MPI
from multiprocessing import Pool


# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# cpus = comm.Get_size()



def log_prior(paras):
    As, omega_m0, omega_bm0_ratio, h = paras
    if 1 < As < 5 and 0.05 < omega_m0 < 0.5 and 0.05 < omega_bm0_ratio < 0.5 and 0.4 < h < 1:
        return 0.0
    else:
        return -numpy.inf


def log_prob(paras, theta_radian, xi, xi_scale, cov_inv, zpts, inv_scale_factor_sq, zhist, z4pk_interp):

    lp = log_prior(paras)

    if not numpy.isfinite(lp):
        # print(lp)
        return -numpy.inf
    else:
        # print(lp)
        # t1 = time.time()
        As, omega_m0, omega_bm0_ratio, h = paras
        omega_bm0 = omega_m0*omega_bm0_ratio
        omega_cm0 = omega_m0 - omega_bm0

        As = As/10**9
        # print(paras)
        try:
            xi_theoretical, sigma8 = cf_tool.get_pk(As, omega_cm0, omega_bm0, h,
                                                      zpts, inv_scale_factor_sq, zhist, z4pk_interp, theta_radian)[:2]

            diff = xi - xi_theoretical.flatten()*xi_scale
            chi_sq = lp - 0.5 * numpy.dot(diff, numpy.dot(cov_inv, diff))
            return chi_sq
        except:
            print(paras)
            return -numpy.inf

        # sigma8_buffer[step_count[0]] = sigma8
        # step_count[0] += 1
        # t2 = time.time()
        # print("%.2f"%(t2-t1), paras,chi_sq)





start = time.time()

expo = int(argv[1])
seed_ini = int(argv[2])
nsteps = int(argv[3])
thread = int(argv[4])
################### read the z data #############################
h5f = h5py.File("./stack_data.hdf5","r")
data = h5f["/data"][()]
h5f.close()

redshift_bin_num = 6
redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)

nz_bin_num = 340

zehist, zebin, zebin_cent = cf_tool.get_nz(data[:, 8], redshift_bin, data[:, 9], nz_bin_num, 2.8)

tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)

inv_scale_factor_sq = (1+zebin_cent)**2

zmin, interp_zmax = 0, 3

z4pk_interp = numpy.linspace(interp_zmax, zmin, 100)



################# read the result data ############################

xi_scale = 1#10**6

# expo_type = "diff_expo"
expo_type = ["diff_expo","same_expo"][expo]


resample_num = 200

npz = numpy.load("./result_cache_%d_%s.npz"%(resample_num,expo_type))
# arcmin to radian
theta_radian = npz["arr_0"]/60/180*numpy.pi
data_num = theta_radian.shape[1]

theta_radian = theta_radian.reshape((tomo_panel_num, int(data_num*1.0/tomo_panel_num)))
xi_p = npz["arr_1"][0]*xi_scale
cov_p = npz["arr_3"]*xi_scale*xi_scale
cov_inv = numpy.linalg.pinv(cov_p)

# prob_coeff = numpy.log(1./(2*numpy.pi)**(data_num/2)/numpy.linalg.det(cov_p)**(0.5))

print("Data vector len: ", theta_radian.shape)



# exit()

################### initialize emcee #############################
numpy.random.seed(seed_ini )#+ rank*10)
nwalkers, ndim = 32, 4
initial = numpy.zeros((nwalkers, ndim))
para_lim = [[1,5],[0.1,0.5],[0.05,0.5],[0.1,1]]
for i in range(ndim):
    a,b = para_lim[i]
    initial[:,i] = numpy.random.uniform(a,b,nwalkers)


# step_count = numpy.zeros((1,),dtype=numpy.intc)
# sigma8_buffer = numpy.zeros((nsteps,))

with Pool(thread) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool, args=(theta_radian, xi_p, xi_scale,
                                                                    cov_inv, zebin_cent,
                                                                    inv_scale_factor_sq, zehist, z4pk_interp))

    sampler.run_mcmc(initial, nsteps, progress=True)

# sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(theta_radian, xi_p,xi_scale,
#                                                                cov_inv, zebin_cent,
#                                                                inv_scale_factor_sq, zehist, z4pk_interp))
# if rank == 0:
#     sampler.run_mcmc(initial, nsteps, progress=True)
# else:
#     sampler.run_mcmc(initial, nsteps, progress=False)
#
# flat_samples = sampler.get_chain(discard=10, thin=1, flat=True)
#
# numpy.save("./data/%s/samples_%d.npz"%(expo_type,rank), flat_samples)

end = time.time()
multi_time = end - start
print("Multiprocessing took {0:.1f} seconds".format(multi_time))
