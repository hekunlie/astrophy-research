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

zmin, interp_zmax = 0, 4

z4pk_interp = numpy.linspace(interp_zmax, zmin, 150)





def log_prob(paras, theta, xip, cov):
    As, omega_m0 = paras

    t = time.time() + numpy.random.uniform(0.05, 0.08)
    while True:
        if time.time() >= t:
            break
    return -0.5 * numpy.sum(theta ** 2)

numpy.random.seed(42)
initial = numpy.random.randn(32, 5)
nwalkers, ndim = initial.shape
nsteps = 1000


with Pool(5) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, log_prob, pool=pool)
    start = time.time()
    sampler.run_mcmc(initial, nsteps, progress=True)
    end = time.time()
    multi_time = end - start
    print("Multiprocessing took {0:.1f} seconds".format(multi_time))
