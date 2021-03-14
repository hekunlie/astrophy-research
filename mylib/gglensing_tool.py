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


def cal_shear(ra, dec, nfw, source_z):
    num = len(ra)
    shear_data = numpy.zeros((4, num))
    for i in range(num):
        kappa = nfw.getConvergence((ra[i], dec[i]), source_z[i])
        g1, g2 = nfw.getLensing((ra[i], dec[i]), source_z[i])[:2]
        shear_data[0, i] = kappa
        shear_data[1, i] = g1
        shear_data[2, i] = g2
        shear_data[3, i] = numpy.arctan2(g2, g1) / 2
    return shear_data
