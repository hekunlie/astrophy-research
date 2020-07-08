import numpy
from scipy import optimize, signal




def convolve(img1, sigma=2):
    my, mx = numpy.mgrid[-2:3, -2:3]
    img2 = numpy.exp(-(my**2+mx**2)/2/sigma**2)/numpy.sqrt(2*numpy.pi)/sigma
    return signal.convolve(img1, img2, mode="same")


def dipole_fun(sin_cos, a1, a2):
    # sin_theta, cos_theta
    return a1 * sin_cos[0] + a2 * sin_cos[1]


def quadrupole_fun(sin_cos_2theta, a1, a2):
    # sin_2theta, cos_2theta
    return a1 * sin_cos_2theta[0] + a2 * sin_cos_2theta[1]

def get_bin_1d(data, bin_num, mode=0, percen_low=0.5, percen_up=99.5, scale=1):

    if mode == 0:
        a = max(numpy.abs(numpy.percentile(data, percen_low)), numpy.percentile(data, percen_up))
        bins = numpy.linspace(-a*scale, a*scale, bin_num + 1)
    #    num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [xy_bin, xy_bin])
    #         print(mode, xbins-ybins)
    else:
        a = numpy.max(numpy.abs(data))
        bins = numpy.linspace(-a*scale, a*scale, bin_num + 1)
    #    num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [xy_bin, xy_bin])
    #         print(mode, xbins-ybins)
    return bins

def get_bin_2d(data_1, data_2, bin_num, mode=0, percen_low=0.5, percen_up=99.5, scale=1):

    if mode == 0:
        a = max(numpy.abs(numpy.percentile(data_1, percen_low)), numpy.percentile(data_1, percen_up))
        b = max(numpy.abs(numpy.percentile(data_2, percen_low)), numpy.percentile(data_2, percen_up))
        xy_max = max(a, b)
        xy_bin = numpy.linspace(-xy_max*scale, xy_max*scale, bin_num + 1)
    #    num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [xy_bin, xy_bin])
    #         print(mode, xbins-ybins)
    else:
        xy_max = max(numpy.max(numpy.abs(data_1)), numpy.max(numpy.abs(data_2)))
        xy_bin = numpy.linspace(-xy_max*scale, xy_max*scale, bin_num + 1)
    #    num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [xy_bin, xy_bin])
    #         print(mode, xbins-ybins)
    return xy_bin


def get_hist_1d(data, bins):
    bin_num = len(bins)-1

    bins_mid = numpy.zeros((bin_num, bin_num))
    num, d1_bins = numpy.histogram(data, bins)

    for i in range(bin_num):
        bins_mid[i] = (bins[i] + bins[i + 1]) / 2

    return num, bins_mid, numpy.abs(bins_mid)


def get_hist_2d(data_1, data_2, bins):
    bin_num = len(bins)-1

    xgrid = numpy.zeros((bin_num, bin_num))
    ygrid = numpy.zeros((bin_num, bin_num))

    num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [bins, bins])

    for i in range(bin_num):
        for j in range(bin_num):
            xgrid[i, j] = (bins[j] + bins[j + 1]) / 2
            ygrid[i, j] = (bins[i] + bins[i + 1]) / 2
    return num, xgrid, ygrid, numpy.sqrt(xgrid ** 2 + ygrid ** 2)


def get_dipole_1d(num, radius, radius_bin_num):

    radius_bin = numpy.linspace(0, radius.max(), radius_bin_num + 1)

    num_dipole = numpy.zeros_like(num)
    raidus_mask = numpy.zeros_like(num)
    mean_num = numpy.zeros_like(num)

    for i in range(radius_bin_num):
        idx1 = radius >= radius_bin[i]
        idx2 = radius < radius_bin[i + 1]
        idx = idx1 & idx2
        mean_of_annuli = num[idx].mean()
        num_dipole[idx] = num[idx] - mean_of_annuli
        mean_num[idx] = mean_of_annuli
        raidus_mask[idx] = i

    return numpy.nan_to_num(num_dipole), radius_bin, raidus_mask, mean_num


def get_dipole_2d(num, radius, radius_bin_num):

    radius_bin = numpy.linspace(0, radius.max(), radius_bin_num + 1)

    num_dipole = numpy.zeros_like(num)
    raidus_mask = numpy.zeros_like(num)
    mean_num = numpy.zeros_like(num)

    for i in range(radius_bin_num):
        idx1 = radius >= radius_bin[i]
        idx2 = radius < radius_bin[i + 1]
        idx = idx1 & idx2
        mean_of_annuli = num[idx].mean()
        num_dipole[idx] = num[idx] - mean_of_annuli
        mean_num[idx] = mean_of_annuli
        raidus_mask[idx] = i

    return numpy.nan_to_num(num_dipole), radius_bin, raidus_mask, mean_num


def get_quadrupole(num_dipole, xgrid, ygrid, radius_bin, radius_bin_num):
    radius_grid = numpy.sqrt(xgrid**2 + ygrid**2)

    sin_theta = ygrid / radius_grid
    cos_theta = xgrid / radius_grid

    num_dipole_fit = numpy.zeros_like(num_dipole)
    num_dipole_fit[:, :] = numpy.nan
    #     print(num_dipole_fit)

    for i in range(radius_bin_num):
        idx1 = radius_grid >= radius_bin[i]
        idx2 = radius_grid < radius_bin[i + 1]
        idx = idx1 & idx2
        if idx.sum() > 5:
            sin_cos = numpy.zeros((2, idx.sum()))
            sin_cos[0] = sin_theta[idx]
            sin_cos[1] = cos_theta[idx]
            res = optimize.curve_fit(dipole_fun, sin_cos, num_dipole[idx])[0]
            num_dipole_fit[idx] = dipole_fun(sin_cos, res[0], res[1])

    num_dipole_fit = numpy.nan_to_num(num_dipole_fit)
    return num_dipole - num_dipole_fit, num_dipole_fit, sin_theta, cos_theta


def fit_quadrupole(quadrupole, xgrid, ygrid, radius_bin, radius_bin_num):
    radius_grid = numpy.sqrt(xgrid ** 2 + ygrid ** 2)
    sin_theta = ygrid / radius_grid
    cos_theta = xgrid / radius_grid
    cos_2theta = cos_theta**2 - sin_theta**2
    sin_2theta = cos_theta*sin_theta*2

    num_quadrupole_fit = numpy.zeros_like(quadrupole)
    num_quadrupole_fit[:, :] = numpy.nan

    for i in range(radius_bin_num):
        idx1 = radius_grid >= radius_bin[i]
        idx2 = radius_grid < radius_bin[i + 1]
        idx = idx1 & idx2
        if idx.sum() > 5:
            sin_cos2 = numpy.zeros((2, idx.sum()))
            sin_cos2[0] = sin_2theta[idx]
            sin_cos2[1] = cos_2theta[idx]
            try:
                res = optimize.curve_fit(quadrupole_fun, sin_cos2, quadrupole[idx])[0]
                num_quadrupole_fit[idx] = quadrupole_fun(sin_cos2, res[0], res[1])
            except:
                print("Failed-fitting")
    return numpy.nan_to_num(num_quadrupole_fit), sin_2theta, cos_2theta
