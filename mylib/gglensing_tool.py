import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
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
        kappa = nfw.getConvergence((ra[i], dec[i]), source_z)
        g1, g2 = nfw.getLensing((ra[i], dec[i]), source_z)[:2]
        shear_data[0, i] = kappa
        shear_data[1, i] = g1
        shear_data[2, i] = g2
        shear_data[3, i] = numpy.arctan2(g2, g1) / 2
    return shear_data


def cal_kappa_map(ra_bin, dec_bin, src_data_list, num_dens, nfw, fq):
    ra_bin_num = len(ra_bin) - 1
    dec_bin_num = len(dec_bin) - 1
    ra_grid, dec_grid = bin2grid(ra_bin, dec_bin)

    kappa_map = numpy.zeros_like(ra_grid)
    kappa_map_err = numpy.zeros_like(ra_grid)

    kappa_model = numpy.zeros_like(ra_grid)

    src_ra, src_dec, src_z, src_mg1, src_mg2, src_mn, src_mu, src_mv = src_data_list

    src_c = SkyCoord(ra=src_ra * units.arcsec, dec=src_dec * units.arcsec, frame="fk5")

    alpha = 1. / num_dens / numpy.pi
    src_num = len(src_ra)

    for i in range(dec_bin_num):
        for j in range(ra_bin_num):
            x, y = ra_grid[i, j], dec_grid[i, j]

            kappa_model[i, j] = nfw.getConvergence((x, y), src_z)

            # separation
            grid_c = SkyCoord(ra=x * units.arcsec, dec=y * units.arcsec, frame="fk5")
            src_sep = grid_c.separation(src_c).arcsec
            src_sep_sq_inv = 1. / src_sep / src_sep
            idx = src_sep_sq_inv >= 0.0001
            pos_ang = grid_c.position_angle(src_c).radian
            print(idx.sum())
            # rotate shear estimators
            mg_t = (src_mg1 * numpy.cos(2 * pos_ang) - src_mg2 * numpy.sin(2 * pos_ang)) * src_sep_sq_inv

            mnu_t = (src_mn - src_mu * numpy.cos(4 * pos_ang) + src_mv * numpy.sin(4 * pos_ang)) * src_sep_sq_inv
            #             print(mg_t.shape, mnu_t.shape)

            k, kerr = fq.find_shear(mg_t[idx], mnu_t[idx], 8)[:2]

            kappa_map[i, j], kappa_map_err[i, j] = numpy.sum(src_sep_sq_inv[idx] * k * alpha), numpy.sum(
                src_sep_sq_inv[idx] * kerr * alpha)

    return kappa_map, kappa_map_err, kappa_model

