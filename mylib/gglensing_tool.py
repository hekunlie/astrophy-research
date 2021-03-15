from sys import path
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py
import galsim
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units


def plot_shear(g1, g2, kappa, ra, dec):
    pass


def bin2grid(xbin, ybin):
    xbin_num = len(xbin) - 1
    ybin_num = len(ybin) - 1
    xgrid = numpy.zeros((ybin_num, xbin_num))
    ygrid = numpy.zeros((ybin_num, xbin_num))

    for i in range(ybin_num):
        for j in range(xbin_num):
            xmid = (xbin[j] + xbin[j + 1]) / 2
            ymid = (ybin[i] + ybin[i + 1]) / 2
            xgrid[i, j] = xmid
            ygrid[i, j] = ymid
    return xgrid, ygrid


def get_pos(ra_bin, dec_bin, bin_num):
    pos_x = numpy.zeros((bin_num, bin_num))
    pos_y = numpy.zeros((bin_num, bin_num))
    for i in range(bin_num):
        for j in range(bin_num):
            x = (ra_bin[i] + ra_bin[i + 1]) / 2
            y = (dec_bin[j] + dec_bin[j + 1]) / 2
            pos_x[i, j] = x
            pos_y[i, j] = y
    return pos_x.reshape((bin_num * bin_num,)), pos_y.reshape((bin_num * bin_num,))


def get_shear(nfw, ra, dec, source_z):
    num = len(ra)
    shear_data = numpy.zeros((5, num))
    for i in range(num):
        kappa = nfw.getConvergence((ra[i], dec[i]), source_z[i])
        gamma1, gamma2 = nfw.getShear((ra[i], dec[i]),source_z[i],reduced=False)
        g1,g2 = gamma1/(1-kappa), gamma2/(1-kappa)
        shear_data[0, i] = kappa
        shear_data[1, i] = gamma1
        shear_data[2, i] = gamma2
        shear_data[3, i] = g1
        shear_data[4, i] = g2
    return shear_data

def get_shear_point(nfw,ra, dec,source_z):
    kappa = nfw.getConvergence((ra, dec), source_z)
    gamma1, gamma2 = nfw.getShear((ra, dec),source_z,reduced=False)
    g1,g2 = gamma1/(1-kappa), gamma2/(1-kappa)
    return kappa, gamma1, gamma2, g1, g2

def get_delta_sigma(nfw, com_dist_len, len_z, com_dist_src, src_z, theta,reduced=False):
    delta_sigma = numpy.zeros_like(theta)
    theta_num = theta.shape[0]
    crit_coeff = 1662895.2081868195*com_dist_src/com_dist_len/(com_dist_src-com_dist_len)/(1+len_z)
    for i in range(theta_num):
        gamma1, gamma2 = nfw.getShear((theta[i],0),src_z,reduced=reduced)
        delta_sigma[i] = numpy.sqrt(gamma1**2 + gamma2**2)*crit_coeff
    return delta_sigma