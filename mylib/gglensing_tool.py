from sys import path
path.append("/home/hklee/work/mylib")
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


def cal_kappa_map(ra_bin, dec_bin, src_data_list, num_dens, nfw):
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


def get_tanshear(data_path, theta_bin_num, bin_min, bin_max, ax, src_num=5000000, msg=False, logbin=False):
    h5f = h5py.File("%s/data_noise_free.hdf5" % data_path, "r")
    data = h5f["/data"][()]
    h5f.close()

    src_all = data.shape[0]
    if src_num >= src_all:
        select_num = src_all
    else:
        select_num = src_num

    mg1 = data[:, 0][:select_num]
    mg2 = data[:, 1][:select_num]
    mn = data[:, 2][:select_num]
    mu = data[:, 3][:select_num]
    mv = data[:, 4][:select_num]

    # shear data
    h5f = h5py.File("%s/source_data.hdf5" % data_path, "r")
    src_ra = h5f["/ra"][()][:select_num]  # arcsec
    src_dec = h5f["/dec"][()][:select_num]  # arcsec
    src_z = h5f["/z"][()][:select_num]
    src_g1 = h5f["/g1"][()][:select_num]
    src_g2 = h5f["/g2"][()][:select_num]
    src_kappa = h5f["/kappa"][()][:select_num]
    src_theta = h5f["/theta"][()][:select_num]
    h5f.close()

    # separation
    len_c = SkyCoord(ra=0 * units.arcsec, dec=0 * units.arcsec, frame="fk5")
    src_c = SkyCoord(ra=src_ra * units.arcsec, dec=src_dec * units.arcsec, frame="fk5")
    src_sep = len_c.separation(src_c).arcsec

    pos_ang = len_c.position_angle(src_c).radian

    # rotate shear estimators
    mg_t = mg1 * numpy.cos(2 * pos_ang) - mg2 * numpy.sin(2 * pos_ang)
    mg_x = mg1 * numpy.sin(2 * pos_ang) + mg2 * numpy.cos(2 * pos_ang)

    mu_t = mu * numpy.cos(4 * pos_ang) - mv * numpy.sin(4 * pos_ang)

    # setup bins
    if logbin:
        theta_bin = tool_box.set_bin_log(bin_min, bin_max, theta_bin_num + 1)
    else:
        theta_bin = numpy.linspace(bin_min, bin_max, theta_bin_num + 1)

    radius_bin = theta_bin / 60 / 60 / 180 * numpy.pi * com_dist_len
    print("Theta_bin: ", theta_bin, " arcsec")
    print("Radius_bin: ", radius_bin, " Mpc/h")
    print("Source plane at Z=%.3f. Source num: %d" % (src_z, select_num))

    theta_pts = numpy.linspace(bin_min, bin_max, 100)
    shear_tan_pts = numpy.zeros_like(theta_pts)
    for i in range(100):
        g1, g2 = nfw.getLensing((theta_pts[i], 0), src_z)[:2]
        shear_tan_pts[i] = numpy.sqrt(g1 ** 2 + g2 ** 2)

    shear_tan = numpy.zeros((5, theta_bin_num))

    for i in range(theta_bin_num):
        idx1 = src_sep >= theta_bin[i]
        idx2 = src_sep < theta_bin[i + 1]
        idx = idx1 & idx2

        shear_tan[0, i] = src_sep[idx].mean()
        if msg:
            print("Theta: [%.5f %.5f], mean: %.5f. Num: %d" % (
            theta_bin[i], theta_bin[i + 1], shear_tan[0, i], idx.sum()))

        mgt_sub = mg_t[idx]
        mgx_sub = mg_x[idx]
        mnut_sub = mn[idx] + mu_t[idx]
        mnux_sub = mn[idx] - mu_t[idx]

        shear_tan[1, i], shear_tan[2, i] = fq.find_shear(mgt_sub, mnut_sub, 8)[:2]
        shear_tan[3, i], shear_tan[4, i] = fq.find_shear(mgx_sub, mnux_sub, 8)[:2]

    ax.plot(theta_pts, shear_tan_pts, label="Model")
    ax.errorbar(shear_tan[0], shear_tan[1], shear_tan[2], label="g_t", mfc="none", marker="s", capsize=4)
    ax.errorbar(shear_tan[0], shear_tan[3], shear_tan[4], label="g_x", mfc="none", marker="s", capsize=4)
    ax.legend()
    xs = ax.set_xlim()
    ax.plot(xs, [0, 0], ls="--", c="gray", alpha=0.5)

    return [src_ra, src_dec, src_z[0], mg1, mg2, mn, mu, mv]