import os
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
from astropy.coordinates import SkyCoord
from astropy import units
from astropy import constants as const


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
separation_bin_num = 3
Rmin, Rmax = 0.05, 0.15 # Mpc/h

separation_bin = hk_tool_box.set_bin_log(Rmin, Rmax, separation_bin_num+1)

# read the parameters
h5f = h5py.File(param_path + "/stack_sheared_para.hdf5", "r")
src_z = h5f["/z"][()]
src_z_m = h5f["/z_m"][()]
src_ra = h5f["/ra"][()]
src_dec = h5f["/dec"][()]
src_g1 = h5f["/gamma1"][()]
src_g2 = h5f["/gamma2"][()]
h5f.close()

gt = numpy.sqrt(src_g1**2 + src_g2**2)
# # the measured ellipticity
src_num = src_g1.shape[0]
rng = numpy.random.RandomState(numpy.random.randint(1, 10000, 1))
e = rng.normal(0, 0.1, src_num)
theta = rng.uniform(0, 2*numpy.pi, src_num)

e1 = numpy.cos(2*theta)
e2 = numpy.sin(2*theta)

src_e1 = e1 + src_g1
src_e2 = e2 + src_g2


# position and separation angle
pos_len = SkyCoord(ra=0 * units.deg, dec=0 * units.deg, frame="fk5")
pos_src = SkyCoord(ra=src_ra * units.deg, dec=src_dec * units.deg, frame="fk5")

separation_radian = pos_len.separation(pos_src).radian
separation_radius = separation_radian * com_dist_len

print("Separation: ",separation_radius.min(), separation_radius.max(),src_ra.max())

position_angle = pos_len.position_angle(pos_src).radian

sin_2theta = numpy.sin(2 * position_angle)
cos_2theta = numpy.cos(2 * position_angle)

# sin_4theta = numpy.sin(4 * position_angle)
# cos_4theta = numpy.cos(4 * position_angle)

src_gt = src_g1 * cos_2theta - src_g2 * sin_2theta
src_gx = src_g1 * sin_2theta + src_g2 * cos_2theta

src_et = src_e1 * cos_2theta - src_e2 * sin_2theta
src_ex = src_e1 * sin_2theta + src_e2 * cos_2theta


h5f = h5py.File(param_path + "/segment_sheared_para.hdf5", "w")
for i in range(separation_bin_num):
    idx1 = separation_radius >= separation_bin[i]
    idx2 = separation_radius < separation_bin[i+1]
    idx = idx1 & idx2
    print("%.4f ~ %.4f Mpc/h %d"%(separation_bin[i], separation_bin[i+1], idx.sum()))

    h5f["/%d/z"%i] = src_z[idx]
    h5f["/%d/z_m"%i] = src_z_m[idx]
    h5f["/%d/ra"%i] = src_ra[idx]
    h5f["/%d/dec"%i] = src_dec[idx]
    h5f["/%d/gamma1"%i] = src_g1[idx]
    h5f["/%d/gamma2"%i] = src_g2[idx]
    h5f["/%d/gamma_t"%i] = src_gt[idx]
    h5f["/%d/gamma_x"%i] = src_gx[idx]
    h5f["/%d/e_t"%i] = src_et[idx]
    h5f["/%d/e_x"%i] = src_ex[idx]
    h5f["/%d/radius"%i] = separation_radius[idx]
    h5f["/%d/radian"%i] = separation_radian[idx]
h5f.close()