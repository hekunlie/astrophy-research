import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import correlation_function_tool as cf_tool
import numpy
import emcee
import time
import h5py
from multiprocessing import Pool



def log_prior(paras):
    As, omega_m0, omega_b_ration = paras
    if 0.01 < As < 5 and 0.05 < omega_m0 < 0.7 and 0.01 < omega_b_ration < 0.5:
        return 0.0
    else:
        return -numpy.inf


def log_prob(paras, theta_radian, xi, cov_inv, zpts, inv_scale_factor_sq, zhist, z4pk_interp):

    lp = log_prior(paras)

    if not numpy.isfinite(lp):
        # print(lp)
        return -numpy.inf
    else:
        # print(lp)
        # t1 = time.time()
        As, omega_m0, omega_b_ration = paras
        omega_bm0 = omega_m0*omega_b_ration
        omega_cm0 = omega_m0 - omega_bm0

        # omega_bm0 = omega_m0*0.02233/0.14213
        # omega_cm0 = omega_m0*0.1198/0.14213
        h = 0.6737
        As = As*10**(-9)
        # print(paras)
        try:
            xi_theoretical, sigma8 = cf_tool.get_tomo_xi(As, omega_cm0, omega_bm0, h,
                                                      zpts, inv_scale_factor_sq, zhist, z4pk_interp, theta_radian)[:2]

            diff = xi - xi_theoretical.flatten()
            chi_sq = lp - 0.5 * numpy.dot(diff, numpy.dot(cov_inv, diff))
            return chi_sq
        except:
            print(As,omega_cm0,omega_bm0)
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
h5f = h5py.File("./data/stack_data.hdf5","r")
data = h5f["/data"][()]
h5f.close()

redshift_bin_num = 5
redshift_bin = numpy.array([0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)

nz_bin_num = 340

zehist, zebin, zebin_cent = cf_tool.get_nz(data[:, 8], redshift_bin, data[:, 9], nz_bin_num, 3)

tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)

inv_scale_factor_sq = (1+zebin_cent)**2

zmin, interp_zmax = 0, 3

z4pk_interp = numpy.linspace(interp_zmax, zmin, 100)


################# read the result data ############################

# expo_type = "diff_expo"
expo_type = ["diff_expo","same_expo"][expo]

resample_num = 200

npz = numpy.load("./data/result_cache_%d_%s.npz"%(resample_num,expo_type))

theta = npz["arr_0"]
# only \xi_+
xi = npz["arr_1"]
cov_inv = npz["arr_4"]

# arcmin to radian
theta_radian = theta/60/180*numpy.pi

data_num = theta_radian.shape[0]
theta_radian = theta_radian.reshape((tomo_panel_num, int(data_num*1.0/tomo_panel_num)))

print("Data vector len: ", xi.shape)


################### initialize emcee #############################
numpy.random.seed(seed_ini)#+ rank*10)
para_lim = [[0.01, 5],[0.01, 0.7],[0.01,0.5]]
nwalkers, ndim = thread, len(para_lim)
initial = numpy.zeros((nwalkers, ndim))

for i in range(ndim):
    a,b = para_lim[i]
    initial[:,i] = numpy.random.uniform(a,b,nwalkers)

print(nsteps, " steps")
with Pool(thread) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool, args=(theta_radian, xi, cov_inv, zebin_cent,
                                                                    inv_scale_factor_sq, zehist, z4pk_interp))

    sampler.run_mcmc(initial, nsteps, progress=True)

# emcee.backends.HDFBackend("./chain_cache.hdf5")

chain = sampler.get_chain()
flat_chain = sampler.get_chain(discard=2000, flat=True)
flat_chain_thin_1 = sampler.get_chain(discard=2000, thin=150, flat=True)

numpy.savez("./data/chain_%s_%d_steps.npz"%(expo_type, nsteps), chain, flat_chain, flat_chain_thin_1)

tau = sampler.get_autocorr_time()
discard_step = int(tau.mean()*2)
thin_step = int(tau.mean()/2)
print(tau, discard_step, thin_step)

flat_chain_thin_2 = sampler.get_chain(discard=discard_step, thin=thin_step, flat=True)

numpy.savez("./data/chain_%s_autocorr_thin_%d_steps.npz"%(expo_type, nsteps), flat_chain_thin_2)

end = time.time()
multi_time = end - start
print("Multiprocessing took {0:.1f} seconds".format(multi_time))
