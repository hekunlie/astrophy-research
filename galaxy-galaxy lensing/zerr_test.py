import os
import matplotlib
from sys import path
path.append("/home/hklee/work/mylib")
from hk_plot_tool import Image_Plot
import hk_tool_box
import hk_gglensing_tool
import numpy
import h5py
import hk_FQlib
import time
# import c4py
import galsim
from astropy.cosmology import FlatLambdaCDM
from mpi4py import MPI



# D(z) = 1/(1+z)*\int_0^z dz/\sqrt{ (1+z)^3\Omega_m + \Omega_Lam}
def integrate(fx, x):
    dx = x[1:] - x[:-1]
    fx_s = (fx[:-1] + fx[1:]) / 2
    return numpy.sum(dx * fx_s)


def get_f_tilde(sigma_z, int_z, inv_sigc, zm):
    sig = sigma_z * (1 + zm)
    fx = numpy.exp(-(int_z - zm) ** 2 / 2 / sig ** 2) / numpy.sqrt(2 * numpy.pi) / sig
    ft = integrate(fx * inv_sigc, int_z)
    dz = numpy.abs(int_z - zm)
    # idx = dz == dz.min()
    #
    # #     print(ft, inv_sigc[idx])
    # #     img = Image_Plot()
    # #     img.subplots(1,1)
    # #     img.axs[0][0].plot(int_z, fx)
    # #     img.show_img()
    # #     img.close_img()
    return ft


def get_bias(sigma_z, dz, pdf_zm, zm, len_z, CF):
    numer = numpy.zeros_like(pdf_zm)
    com_dist_src = CF.com_distance(zm)

    zt = numpy.linspace(0.001, 5, 501)
    com_dist_t = CF.com_distance(zt)

    inv_sigc = 1. / CF.get_sigma_crit(zt)
    idx = zt <= len_z
    inv_sigc[idx] = 0

    for i in range(len(zm)):
        numer[i] = get_f_tilde(sigma_z, zt, inv_sigc, zm[i]) * pdf_zm[i]

    inv_sigc_zm = 1. / CF.get_sigma_crit(zm)
    idx = zm <= len_z
    inv_sigc_zm[idx] = 0
    denorm = pdf_zm * inv_sigc_zm

    idx = zm <= len_z + dz
    numer[idx] = 0
    denorm[idx] = 0
    numer_int = integrate(numer, zm)
    denorm_int = integrate(denorm, zm)
    #     print(denorm_int)

    bias = numer_int / denorm_int

    #     print(numer)
    #     print(denorm)
    return bias


def get_inv_sigma_2prime(len_z, com_dist_len, src_z, com_dist_src, omega_m):
    beta = (1 + len_z) * com_dist_len ** 2 / 1662916.5401756007
    c_over_H = 2997.92458
    fz = (1 + src_z) ** 3 * omega_m + 1 - omega_m
    inv_sigc_2p = -beta * c_over_H * ( 2 / com_dist_src ** 3 * c_over_H / fz + 1.5 * (1 + src_z) ** 2 * omega_m / com_dist_src ** 2 / fz ** (
            1.5))
    idx = src_z < len_z
    inv_sigc_2p[idx] = 0
    return inv_sigc_2p


def get_inv_sigma_prime(len_z, com_dist_len, src_z, com_dist_src, omega_m, order):
    beta = (1 + len_z) * com_dist_len ** 2 / 1662916.5401756007
    c_over_H = 2997.92458
    fz = (1 + src_z) ** 3 * omega_m + 1 - omega_m

    w_p = c_over_H * fz ** (-0.5)
    w_2p = -c_over_H * 1.5 * omega_m * (1 + src_z) ** 2 * fz ** (-1.5)
    w_3p = -c_over_H * 1.5 * omega_m * (3 * omega_m * (1 + src_z) ** 4 * fz ** (-2.5) + 2 * (1 + src_z) * fz ** (-1.5))

    if order == 1:
        inv_sigc_p = beta / com_dist_src ** 2 * w_p
    elif order == 2:
        inv_sigc_p = beta * (-2 * w_2p ** 2 / com_dist_src ** 3 + w_2p / com_dist_src ** 2)
    elif order == 3:
        inv_sigc_p = beta * (6 * w_p ** 2 / com_dist_src ** 4 - 4 * w_p * w_2p / com_dist_src ** 3 - 2 * w_2p / com_dist_src * 3 + w_3p / com_dist_src ** 2)

    idx = src_z <= len_z
    inv_sigc_p[idx] = 0
    return inv_sigc_p


def get_bias_approx(sigma_z, dz, pdf_zm, zm, len_z, CF):
    com_dist_len = CF.com_distance(len_z)
    com_dist_src = CF.com_distance(zm)

    inv_sigc = 1. / CF.get_sigma_crit(zm)
    inv_sigc_2p = get_inv_sigma_2prime(len_z, com_dist_len, zm, com_dist_src, CF.Omega_m0)

    numer = pdf_zm * inv_sigc_2p
    denorm = pdf_zm / CF.get_sigma_crit(zm)

    idx = zm <= len_z + dz
    numer[idx] = 0
    denorm[idx] = 0
    numer_int = (numer.sum() - numer[0] / 2 - numer[-1] / 2) * (zm[1] - zm[0])
    denorm_int = (denorm.sum() - denorm[0] / 2 - denorm[-1] / 2) * (zm[1] - zm[0])
    #     print(denorm_int)
    bias = 1 + 0.5 * sigma_z ** 2 * numer_int / denorm_int
    return bias


def get_bias_approx_new(sigma_z, dz, pdf_zm, zm, len_z, CF):
    numer = numpy.zeros_like(pdf_zm)

    zt = numpy.linspace(0.001, 5, 501)

    inv_sigc_zm = 1. / CF.get_sigma_crit(zm)
    idx = zm <= len_z
    inv_sigc_zm[idx] = 0

    for i in range(len(zm)):
        numer[i] = get_f_tilde(sigma_z, zt, inv_sigc_zm, zm[i]) * pdf_zm[i]

    denorm = pdf_zm * inv_sigc_zm

    idx = zm <= len_z + dz
    numer[idx] = 0
    denorm[idx] = 0
    numer_int = integrate(numer, zm)
    denorm_int = integrate(denorm, zm)
    bias = numer_int / denorm_int
    return bias


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


# cosmology
omega_m0 = 0.31
omega_lam0 = 1 - omega_m0
h = 0.6735
H_0 = 100 * h

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)

# Halo parameters
Mass = 3*10 ** 13  # M_sun/h
conc = 6  # concentration
len_z = 0.2  # redshift
halo_position = galsim.PositionD(0, 0)  # arcsec
com_dist_len = cosmos.comoving_distance(len_z).value * h  # Mpc/h
print("Lens plane at z = %.2f, %.5f Mpc/h" % (len_z, com_dist_len))

# lens profile
CF = hk_gglensing_tool.Cosmos_flat(omega_m0, 100*h)
CF.NFW_profile_galsim((0,0), Mass, conc, len_z)

param_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_9/params"
separation_bin_num = 1
Rmin, Rmax = 0.05, 0.07  # Mpc/h

separation_bin = hk_tool_box.set_bin_log(Rmin, Rmax, separation_bin_num + 1)

result = numpy.zeros((6, separation_bin_num))

zpts = numpy.linspace(0, 2, 301)
zpts_mid = (zpts[1:] + zpts[:-1]) / 2
# print(zpts)
dz = 0.2
sigma_z = 0.05
files = ["segment_sheared_para_mean_0.4_errsig_z0.05.hdf5",
         "segment_sheared_para_mean_0.6_errsig_z0.05.hdf5",
         "segment_sheared_para_mean_0.8_errsig_z0.05.hdf5"]

ii = rank
radius_tag = 0
# read the parameters
h5f = h5py.File(param_path + "/%s" % files[ii], "r")
src_z = h5f["/%d/z" % radius_tag][()]
src_z_m = h5f["/%d/z_m" % radius_tag][()]
src_radius = h5f["/%d/radius" % radius_tag][()]
src_radian = h5f["/%d/radian" % radius_tag][()]
src_et = h5f["/%d/e_t" % radius_tag][()]
h5f.close()

pz = numpy.histogram(src_z, zpts)[0]
pz_norm = pz / pz.sum() / (zpts[1] - zpts[0])
pz_m = numpy.histogram(src_z_m, zpts)[0]
pz_m_norm = pz_m / pz_m.sum() / (zpts[1] - zpts[0])

# if radius_tag == 0:
#     img = Image_Plot()
#     img.subplots(1, 3)
#     img.axs[0][0].hist(src_z, 100, histtype="step", label="$z$")
#     img.axs[0][0].hist(src_z_m, 100, histtype="step", label="$z_m$")
#     idx1 = src_z_m > 0.5
#     idx2 = src_z_m < 0.51
#     idx = idx1 & idx2
#     img.axs[0][1].hist(src_z[idx], 100, histtype="step", label="$z_m$")
#     img.axs[0][2].plot(zpts_mid, pz)
#     img.axs[0][2].plot(zpts_mid, pz_m)
#
#     img.axs[0][0].legend()
#     img.show_img()

idx_z = src_z > len_z + dz
idx_zm = src_z_m > len_z + dz

bias_approx = get_bias_approx(sigma_z, dz, pz_m_norm, zpts_mid, len_z, CF)
bias = get_bias(sigma_z, dz, pz_m_norm, zpts_mid, len_z, CF)

src_num = src_et.shape[0]

sigma_crit = CF.get_sigma_crit(src_z)
src_ds = src_et * sigma_crit

sigma_crit_m = CF.get_sigma_crit(src_z_m)
src_ds_m = src_et * sigma_crit_m

sep_arcsec = src_radian[idx_z].mean() / numpy.pi * 180 * 3600
sep_arcsec_m = src_radian[idx_zm].mean() / numpy.pi * 180 * 3600

ds_true = CF.get_shear(sep_arcsec, 0, src_z[idx_z].mean())[2]
img = Image_Plot()
img.subplots(1, 1)

#     ds1, ds1_err = hk_FQlib.find_shear_cpp(src_ds_m[idx_zm], numpy.ones_like(src_ds_m)[idx_zm],
#                                                                    numpy.ones_like(src_ds_m)[idx_zm], 10,
#                                                                    scale =1, left=-100, right=300,chi_gap=20., max_iters=50)[:2]
#     print("True: %.4f. Measured: %.4f(%.4f). Ratio: %.4f. Expected Bias: %.4f %.4f"%(ds_true, ds1, ds1_err, ds1/ds_true, bias, bias_approx))
ds2, ds2_err = hk_FQlib.find_shear_cpp(src_et[idx_zm], numpy.ones_like(src_et[idx_zm]) / sigma_crit_m[idx_zm],
                                       numpy.ones_like(src_et[idx_zm]), 8,
                                       fit_scale=1., left=-200, right=500, chi_gap=20., max_iters=50,
                                       fig_ax=img.axs[0][0])[:2]
print("True: %.4f. Measured: %.4f(%.4f). Ratio: %.4f(%.4f). Expected Bias: %.4f %.4f" % (
ds_true, ds2, ds2_err, ds2 / ds_true, ds2_err / ds_true, bias, bias_approx))
img.show_img()
h5f = h5py.File(param_path + "/ratio_z_%s" % files[ii], "w")
h5f["/%d/result" % radius_tag] = numpy.array([ds_true, ds2, ds2_err, ds2 / ds_true, ds2_err / ds_true, bias, bias_approx])
h5f.close()
