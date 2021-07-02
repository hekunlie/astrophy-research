import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import correlation_function_tool as cf_tool
import tool_box
import numpy
import emcee
import time
import h5py
from multiprocessing import Pool


def log_prior_ccl(paras):
    sigma8, omega_cm0, omega_bm0 = paras
    if 0.1 < sigma8 < 3 and 0.01 < omega_cm0 < 0.8 and 0.01 < omega_bm0 < 0.2:
        return 0.0
    else:
        return -numpy.inf


def log_prob_ccl(paras):
    # theta_deg, xi, cov_inv, theta_num_per_bin, zpts, zhist, ell, used_zbin):

    lp = log_prior_ccl(paras)

    if not numpy.isfinite(lp):
        # print(lp)
        return -numpy.inf
    else:
        # print(lp)
        # t1 = time.time()
        sigma8, omega_cm0, omega_bm0 = paras

        h = 0.674

        # print(paras)
        try:
            xi_predict = cf_tool.get_tomo_xi_ccl_2(sigma8, omega_cm0, omega_bm0, h, zpts, zehist,
                                                   theta_degree, theta_num_per_bin, ell, used_zbin)

            diff = xi - xi_predict
            chi_sq = lp - 0.5 * numpy.dot(diff, numpy.dot(cov_inv, diff))

            return chi_sq
        except:
            print(sigma8, omega_cm0, omega_bm0)
            return -numpy.inf

def log_prior_ccl_fixb(paras):
    sigma8, omega_cm0 = paras
    if 0.1 < sigma8 < 3 and 0.01 < omega_cm0 < 1:
        return 0.0
    else:
        return -numpy.inf


def log_prob_ccl_fixb(paras):
    # theta_deg, xi, cov_inv, theta_num_per_bin, zpts, zhist, ell, used_zbin):

    lp = log_prior_ccl_fixb(paras)

    if not numpy.isfinite(lp):
        # print(lp)
        return -numpy.inf
    else:
        # print(lp)
        # t1 = time.time()
        sigma8, omega_cm0 = paras
        omega_bm0 = 0.049
        h = 0.674

        # print(paras)
        try:
            xi_predict = cf_tool.get_tomo_xi_ccl_2(sigma8, omega_cm0, omega_bm0, h, zpts, zehist,
                                                   theta_degree, theta_num_per_bin, ell, used_zbin)

            diff = xi - xi_predict
            chi_sq = lp - 0.5 * numpy.dot(diff, numpy.dot(cov_inv, diff))

            return chi_sq
        except:
            print(sigma8, omega_cm0, omega_bm0)
            return -numpy.inf


start = time.time()

expo = int(argv[1])
seed_ini = int(argv[2])
nsteps = int(argv[3])
thread = int(argv[4])


################### read the z data #############################

redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)

h5f = h5py.File("./data/zhist.hdf5","r")
zehist = h5f["/zhist"][()]
zebin = h5f["/zbin"][()]
zpts = h5f["/zbin_cent"][()]
h5f.close()
redshift_bin_num = 6
tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)


################# read the result data ############################

expo_type = ["diff_expo","same_expo","stack_expo"][expo]

resample_num = 200

h5f = h5py.File("./data/result_cache_%d_%s.hdf5"%(resample_num,expo_type),"r")

theta = h5f["/theta"][()]
# only \xi_+
xi = h5f["/xi_p"][()]
cov_inv = h5f["/inv_cov_xi_p"][()]
used_zbin = h5f["/used_zbin"][()]
h5f.close()

# arcmin to radian & degree
theta_radian = theta/60/180*numpy.pi
theta_degree = theta/60

data_num = theta_radian.shape[0]

theta_num_per_bin = 5

ell = tool_box.set_bin_log(5, 15000, 200)

print("Data vector len: ", xi.shape)

################### initialize emcee #############################
numpy.random.seed(seed_ini)#+ rank*10)
para_lim = [[0.1, 3],[0.01, 1]]#,[0.01,0.2]]
nwalkers, ndim = thread, len(para_lim)
initial = numpy.zeros((nwalkers, ndim))

for i in range(ndim):
    a,b = para_lim[i]
    initial[:,i] = numpy.random.uniform(a,b,nwalkers)

print(nsteps, " steps")
with Pool(thread) as pool:
    # sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_ccl, pool=pool)
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob_ccl_fixb, pool=pool)

    sampler.run_mcmc(initial, nsteps, progress=True)


chain = sampler.get_chain()
flat_chain = sampler.get_chain(discard=1000, flat=True)
flat_chain_thin_1 = sampler.get_chain(discard=1000, thin=10, flat=True)

numpy.savez("./chain/chain_%s_%d_steps.npz"%(expo_type, nsteps), chain, flat_chain, flat_chain_thin_1)

tau = sampler.get_autocorr_time()
discard_step = int(tau.mean()*2)
thin_step = int(tau.mean()/2)
print(tau, discard_step, thin_step)

flat_chain_thin_2 = sampler.get_chain(discard=discard_step, thin=thin_step, flat=True)

numpy.savez("./chain/chain_%s_autocorr_thin_%d_steps.npz"%(expo_type, nsteps), flat_chain_thin_2)

end = time.time()
multi_time = end - start
print("Multiprocessing took {0:.1f} seconds".format(multi_time))
