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

z4pk_interp = numpy.linspace(0, interp_zmax, 150)


def get_pk(As, Omega_cm0, Omega_bm0, h, zpts, inv_scale_factor_sq, zhist, z4pk_interp, theta_num, theta_radian):
    c = 299792.458  # Km/s
    ns = 0.965
    H0 = 100*h
    omg_cm0h2 = Omega_cm0 * h ** 2

    omg_bm0h2 = Omega_bm0 * h ** 2

    omg_m0 = Omega_cm0 + Omega_bm0

    alpha = 9 / 4 * omg_m0 ** 2 * (100 / c) ** 4  # h^4 \cdot Mpc^{-4}

    redshift_bin_num = 6
    zpts_num = 400

    Lpts_min, Lpts_max, Lpts_num = 10, 7000, 6000

    kmin, cal_kmax, interp_kmax, kpts_num = 1e-4, 3, 4, 300

    # tomo panel num
    tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)

    # comoving distance Mpc/h
    com_dist = cf_tool.z2dist(zpts, H0, omg_m0)
    delta_com_dist = com_dist[1:] - com_dist[:-1]

    # lense efficiency
    qx = cf_tool.get_lenq(zhist, com_dist, tomo_panel_num, zpts_num)

    # g(x)^2/a(x)^2 in the integration
    integ_factor = numpy.zeros((tomo_panel_num, zpts_num))
    tag = 0
    for i in range(redshift_bin_num):
        for j in range(i, redshift_bin_num):
            integ_factor[tag] = qx[i] * qx[j] * inv_scale_factor_sq
            tag += 1

    # set up bins for L of P(L), and decide the ks needed
    Lpts, dLpts, integ_kh = cf_tool.ready4PL(Lpts_min, Lpts_max, Lpts_num, kmin, cal_kmax, com_dist, zpts_num)

    Lpts_theta = numpy.zeros((theta_num, Lpts_num))

    for i in range(theta_num):
        Lpts_theta[i] = theta_radian[i] * Lpts

    # Pk interpolation
    integ_pk = numpy.zeros((zpts_num, Lpts_num))
    PLs = numpy.zeros((tomo_panel_num, Lpts_num))
    # the theoretical line
    xi_plus = numpy.zeros((tomo_panel_num, theta_num))

    camb_result = cf_tool.get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, z4pk_interp, kmax=interp_kmax)[0]
    pk_interp = camb_result.get_matter_power_interpolator(nonlinear=True, hubble_units=True, k_hunit=True)
    sigma8_now = camb_result.get_sigma8()[-1]
    print("sigma8: %.4f" % sigma8_now)

    for i in range(zpts_num):
        idx = integ_kh[i] > 0
        integ_pk[i][idx] = pk_interp.P(zpts[i], integ_kh[i][idx])

    # get P(L)
    for i in range(tomo_panel_num):
        pls_i = cf_tool.get_PL(integ_pk, Lpts_num, integ_factor[i], delta_com_dist) * alpha
        PLs[i] = pls_i

        # calculate the \chi_plus(\theta) , the theoretical line
        for j in range(theta_num):
            integ_part = Lpts * scipy.special.j0(Lpts_theta[j]) * PLs[i]
            xi_plus[i, j] = numpy.sum(((integ_part[1:] + integ_part[:-1]) / 2 * dLpts)) / 2 / numpy.pi

    return xi_plus


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
