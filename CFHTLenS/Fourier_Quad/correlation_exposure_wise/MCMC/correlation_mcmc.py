import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import tool_box
import correlation_function_tool as cf_tool
from multiprocessing import Pool
import time
import numpy
import emcee
import time
import scipy
import h5py


def log_prior(paras):
    As, omega_m0, omega_bm0_ratio, h = paras
    if 0 < As < 5 and 0 < omega_m0 < 0.5 and 0 < omega_bm0_ratio < 1 and 0.3 < h < 1:
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
        t1 = time.time()
        As, omega_m0, omega_bm0_ratio, h = paras
        omega_bm0 = omega_m0*omega_bm0_ratio
        omega_cm0 = omega_m0 - omega_bm0

        As = As/10**9
        # print(paras)
        xi_theoretical, sigma8 = cf_tool.get_pk(As, omega_cm0, omega_bm0, h,
                                                  zpts, inv_scale_factor_sq, zhist, z4pk_interp, theta_radian)[:2]

        diff = xi - xi_theoretical*xi_scale
        print(diff.shape,cov_inv.shape)
        # sigma8_buffer[step_count[0]] = sigma8
        # step_count[0] += 1
        t2 = time.time()
        print("%.2f"%(t2-t1), paras)
        return lp - 0.5*numpy.dot(diff, numpy.dot(cov_inv, diff))



################# read the result data ############################

xi_scale = 10**6

expo_type = "diff_expo"
# expo_type = "same_expo"

resample_num = 200

npz = numpy.load("./result_cache_%d_%s.npz"%(resample_num,expo_type))
# arcmin to radian
theta_radian = npz["arr_0"][0]/60/180*numpy.pi
xi_p = npz["arr_1"][0]*xi_scale
cov_p = npz["arr_3"]*xi_scale*xi_scale
cov_inv = numpy.linalg.inv(cov_p)

data_num = theta_radian.shape[0]
# prob_coeff = numpy.log(1./(2*numpy.pi)**(data_num/2)/numpy.linalg.det(cov_p)**(0.5))

print("Data vector len: ", data_num)
print(theta_radian.shape,xi_p.shape)
# print(cov_p)
# print(numpy.linalg.det(cov_p))
# print(numpy.dot(xi_p,numpy.dot(cov_inv, xi_p)))

# exit()
################### read the z data #############################
h5f = h5py.File("/mnt/perc/hklee/CFHT/correlation/cata/stack_data.hdf5","r")
data = h5f["/data"][()]
h5f.close()

redshift_bin_num = 6
redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)

nz_bin_num = 400

zehist, zebin, zebin_cent = cf_tool.get_nz(data[:, 8], redshift_bin, data[:, 9], nz_bin_num, 2.2)

inv_scale_factor_sq = (1+zebin_cent)**2

zmin, interp_zmax = 0, 3

z4pk_interp = numpy.linspace(interp_zmax, zmin, 150)



################### initialize emcee #############################
numpy.random.seed(12312)
nwalkers, ndim = 32,4
initial = numpy.random.randn(nwalkers, ndim)#.reshape((nwalkers, ndim))

nsteps = 1000
step_count = numpy.zeros((1,),dtype=numpy.intc)
sigma8_buffer = numpy.zeros((nsteps,))


sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, args=(theta_radian, xi_p,xi_scale,
                                                               cov_inv, zebin_cent,
                                                               inv_scale_factor_sq, zehist, z4pk_interp
                                                               ))

start = time.time()

sampler.run_mcmc(initial, nsteps, progress=True)

flat_samples = sampler.get_chain(discard=10, thin=15, flat=True)

numpy.save("./samples.npz", flat_samples)

end = time.time()
multi_time = end - start
print("Multiprocessing took {0:.1f} seconds".format(multi_time))
