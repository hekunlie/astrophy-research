from sys import path
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import gglensing_tool
import numpy
import h5py
import galsim
import FQlib
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
import time


def get_diff(hist2d):
    hist_num = hist2d.shape[1]
    hist_num2 = int(hist_num / 2)

    hist2d_left = numpy.flip(hist2d[:, :hist_num2], axis=1)
    hist2d_right = hist2d[:, hist_num2:]

    hist2d_diff = hist2d_right - hist2d_left
    hist2d_sum = hist2d_right + hist2d_left

    idx = hist2d_sum > 0
    chisq = numpy.zeros_like(hist2d_left)
    chisq[idx] = hist2d_diff[idx] ** 2 / hist2d_sum[idx] / 2

    return chisq, hist2d_left, hist2d_right, hist2d_diff, hist2d_sum


def prepare_data(data_path, para_path, cosmos, len_z, H_0, len_pos, z_err):
    h5f = h5py.File(data_path, "r")
    data = h5f["/data"][()]
    h5f.close()

    h5f = h5py.File(para_path, "r")
    ra = h5f["/ra"][()]
    dec = h5f["/dec"][()]
    z = h5f["/z"][()]
    h5f.close()

    com_dist = cosmos.comoving_distance(z+z_err).value * H_0 / 100
    com_dist_len = cosmos.comoving_distance(len_z).value * H_0 / 100

    crit_coeff = 1662895.2081868195 * com_dist / com_dist_len / (com_dist - com_dist_len) / (1 + len_z)

    src_pos = SkyCoord(ra=ra * units.arcsec, dec=dec * units.arcsec, frame="fk5")
    position_theta = len_pos.position_angle(src_pos).rad

    mg1r, mg2r, mnur1, mnur2 = FQlib.rotate(data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], position_theta)
    mg1r *= crit_coeff

    return mg1r, mnur1


# cosmology
omega_m0 = 0.31
omega_lam0 = 1 - omega_m0
h = 0.6735
C_0_hat = 2.99792458
H_0 = 100*h
coeff = 1000*C_0_hat/h

coeff_crit = C_0_hat**2/4/numpy.pi/6.674

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)


# Halo parameters
Mass = 5*10**12.5 #M_sun/h
conc = 3.5 #concentration
len_z = 0.3 #redshift
halo_position = galsim.PositionD(0,0) #arcsec
com_dist_len = cosmos.comoving_distance(len_z).value*h #Mpc/h
print("Lens plane at z = %.2f, %.5f Mpc/h"%(len_z, com_dist_len))
len_pos = SkyCoord(ra=0*units.arcsec, dec=0*units.arcsec,frame="fk5")

# lens profile
nfw = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)


# data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z"
data_path = "/coma/hklee/Galaxy_Galaxy_lensing_test/background/continue_source_z"

hist_bin_tag = 0
shear_tag = 0
foreground_z_err = numpy.abs(numpy.random.normal(0, 0.2, 20000000)) + 0.31
num_total = 20000000
num_s = 16000000
num_non = num_total - num_s

dilution_ratio = num_non*1.0/num_total

data_path1 = data_path + "/data/sheared_data_%d_noise_free.hdf5"%shear_tag
data_path2 = data_path + "/data/non_sheared_data_%d_noise_free.hdf5"%shear_tag
para_path1 = data_path + "/paras/sheared_para_%d.hdf5"%shear_tag
para_path2 = data_path + "/paras/non_sheared_para_%d.hdf5"%shear_tag

mg1r_src, mnur1_src = prepare_data(data_path1, para_path1, cosmos, len_z, H_0, len_pos, 0)

mg1r_non, mnur1_non = prepare_data(data_path2, para_path2, cosmos, len_z, H_0, len_pos, foreground_z_err)


mg1r_total = numpy.zeros((num_total, ))
mnur1_total = numpy.zeros((num_total, ))

mg1r_total[:num_s] = mg1r_src[:num_s]
mg1r_total[num_s:] = mg1r_non[:num_non]

mnur1_total[:num_s] = mnur1_src[:num_s]
mnur1_total[num_s:] = mnur1_non[:num_non]


# set bins for 2d histogram
hist2d_bin_num = 1000
hist2d_bin_num2 = int(hist2d_bin_num/2)
mg_bin = gglensing_tool.set_bin(mg1r_src, hist2d_bin_num, bound_scale=1.1, method="log")
mnur_bin = gglensing_tool.set_bin(mnur1_src, hist2d_bin_num, bound_scale=1.1, method="log")

x = numpy.zeros((hist2d_bin_num,hist2d_bin_num))
y = numpy.zeros((hist2d_bin_num,hist2d_bin_num))

for i in range(hist2d_bin_num):
    x[:,i] = (mg_bin[i] + mg_bin[i+1])/2
    y[i] = (mnur_bin[i] + mnur_bin[i+1])/2

xh, yh = x[:,hist2d_bin_num2:], y[:,hist2d_bin_num2:]


img = Image_Plot(xpad=0.25, ypad=0.25)
img.subplots(1,2)

t1 = time.time()

src_result = FQlib.find_shear_cpp(mg1r_src, mnur1_src, 20, left=-100, right=200,
                                  fit_num=30, chi_gap=50,fig_ax=img.axs[0][0])
src_Ds, src_Ds_err, src_coeffs, src_chisqs_min, src_bins = src_result

t2 = time.time()

print("%d foreground. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_s,src_Ds, src_Ds_err, src_chisqs_min, t2-t1))

non_result = FQlib.find_shear_cpp(mg1r_non, mnur1_non, 20, left=-100, right=200,
                                  fit_num=30, chi_gap=50,fig_ax=img.axs[0][1])
non_Ds, non_Ds_err,non_coeffs, non_chisqs_min, non_bins = non_result

t3 = time.time()
print("%d foreground. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_s, non_Ds, non_Ds_err, non_chisqs_min, t3-t2))
img.save_img("./pdf_sym.png")
# img.show_img()

