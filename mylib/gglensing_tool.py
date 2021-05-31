from sys import path
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py


import ctypes
import numpy.ctypeslib as ctl
histlib = ctypes.cdll.LoadLibrary("/home/hklee/work/mylib/libc4py.so")

hist2d_fast = histlib.hist2d_fast
hist2d_fast.restype = None
hist2d_fast.argtypes = [ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                                 ctypes.c_int,
                                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                                 ctypes.c_int,
                                 ctypes.c_int,
                                 ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous'),
                                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous')]




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


def asymetric(num_hist):
    bin_num = num_hist.shape[0]
    bin_num2 = int(bin_num / 2)
    inverse = range(bin_num2 - 1, -1, -1)

    left_bin = num_hist[:bin_num2][inverse]
    right_bin = num_hist[bin_num2:]

    chisq = (left_bin ** 2 + right_bin ** 2) / 2 / left_bin / right_bin
    chisq = numpy.sum(chisq) - bin_num2

    return chisq


def data_mix(data_path1, para_path1, data_path2, para_path2, num1, num2, dz):
    data_mix = numpy.zeros((num1 + num2, 5))
    ra_mix, dec_mix, z_mix = numpy.zeros((num1 + num2,)), numpy.zeros((num1 + num2,)), numpy.zeros((num1 + num2,))

    h5f = h5py.File(data_path1, "r")
    data1 = h5f["/data"][()]
    h5f.close()
    h5f = h5py.File(data_path2, "r")
    data2 = h5f["/data"][()]
    h5f.close()

    data_mix[:num1] = data1[:num1]
    data_mix[num1:] = data2[:num2]

    h5f = h5py.File(para_path1, "r")
    ra1 = h5f["/ra"][()]
    dec1 = h5f["/dec"][()]
    z1 = h5f["/z"][()]
    h5f.close()
    h5f = h5py.File(para_path2, "r")
    ra2 = h5f["/ra"][()]
    dec2 = h5f["/dec"][()]
    z2 = h5f["/z"][()]
    h5f.close()

    ra_mix[:num1] = ra1[:num1]
    dec_mix[:num1] = dec1[:num1]
    z_mix[:num1] = z1[:num1]

    ra_mix[num1:] = ra2[:num2]
    dec_mix[num1:] = dec2[:num2]
    z_mix[num1:] = z2[:num2] + dz

    return data_mix[:, 0], data_mix[:, 1], data_mix[:, 2], data_mix[:, 3], data_mix[:, 4], ra_mix, dec_mix, z_mix


def show_signal(radius, signals, plot_scale=100):
    img = Image_Plot(xpad=0.25, ypad=0.15)
    img.subplots(1, 2)
    img.axs[0][0].errorbar(radius, signals[0] * plot_scale, signals[1] * plot_scale, capsize=3, marker="s", c="r",
                           fmt="--", label="$\gamma_t$")
    img.axs[0][1].errorbar(radius, signals[2] * plot_scale, signals[3] * plot_scale, capsize=3, marker='s', c="r",
                           fmt="--", label="$\gamma_x$")
    img.axs[0][0].errorbar(radius, true_signal * plot_scale, c="b", fmt="--", label="true $\gamma_t$")

    img.axs[0][0].set_yscale("log")

    for i in range(2):
        img.axs[0][i].legend(fontsize=img.legend_size)
        img.axs[0][i].set_xscale("log")

        img.set_label(0, i, 0, "$10^2\gamma$")
        img.set_label(0, i, 1, "Radius [Mpc/h]")
    ys = img.axs[0][1].set_ylim()
    xs = img.axs[0][1].set_xlim()
    img.axs[0][1].plot([xs[0], xs[1]], [0, 0], ls="--", c="k")
    xs = img.axs[0][1].set_xlim()
    img.axs[0][1].set_ylim(ys)
    img.show_img()


def set_bin(data, bin_num, bound_scale, method="log", log_end=5):
    if method == "mean":
        temp_data = numpy.sort(data)
        bin_size = len(temp_data) / bin_num
        bins = [temp_data[int(i * bin_size)] for i in range(bin_num)]
        bins.append(temp_data[-1] * bound_scale)
        bins = numpy.array(bins)
    elif method == "sym-mean":
        temp_data = numpy.sort(numpy.abs(data))  # [:int(len(data[data>0])*0.99)]
        bin_size = len(temp_data) / bin_num * 2
        bins = numpy.array([temp_data[int(i * bin_size)] for i in range(1, int(bin_num / 2))])
        bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
        bound = numpy.max(numpy.abs(data)) * bound_scale
        bins = numpy.append(-bound, numpy.append(bins, bound))

    elif method == "mix_bin":
        mean_num, log_num = bin_num
        mean_num2, log_num2 = int(mean_num / 2), int(log_num / 2)

        hbins = numpy.zeros((mean_num2 + log_num2,))

        temp_data = numpy.sort(numpy.abs(data))  # [:int(len(data[data>0])*0.99)]

        bin_size = len(temp_data) / mean_num * 2
        mean_bins = numpy.array([temp_data[int(i * bin_size)] for i in range(1, mean_num2)])

        hbins[:mean_num2 - 1] = mean_bins

        hbins[mean_num2 - 1] = temp_data[-1]

        ed = numpy.log10(temp_data[-1])
        log_bin = 10 ** numpy.linspace(ed, ed + log_end, log_num2 + 1)
        hbins[mean_num2:] = log_bin[1:]

        inv = range(mean_num2 + log_num2 - 1, -1, -1)
        bins = numpy.zeros((mean_num + log_num + 1))
        bins[mean_num2 + log_num2 + 1:] = hbins
        bins[:mean_num2 + log_num2] = -hbins[inv]

    elif method == "log":
        data_min = data.min()
        bin_num2 = int(bin_num / 2)
        bins = numpy.zeros((bin_num + 1,))

        if data_min < 0:
            temp_data = numpy.sort(numpy.abs(data))
            data_min, data_max = temp_data[0], temp_data[-1]
            #             print(data_min, data_max)
            if data_min < 0.1:
                data_min = 0.1
            else:
                data_min = data_min * 0.95
            bin_num_ = bin_num2 - 1
            inverse = range(bin_num_, -1, -1)
            hbins = tool_box.set_bin_log(data_min, data_max, bin_num2)

            #             hbins = numpy.exp(numpy.linspace(numpy.log(data_min), numpy.log(data_max), bin_num2))
            #             hbins = 10**numpy.linspace(numpy.log10(data_min), numpy.log10(data_max), bin_num2)

            bins[bin_num2 + 1:] = hbins
            bins[:bin_num2] = -hbins[inverse]

            bins[0] = bins[0] * bound_scale
            bins[-1] = bins[-1] * bound_scale
        else:
            data_max = data.max()
            bins = tool_box.set_bin_log(data_min, data_max, bin_num+1)
            bins[0] = bins[0] * 0.95
            bins[-1] = bins[-1] * bound_scale
    else:
        m = max(numpy.abs(data.min()), numpy.abs(data.max()))
        bins = numpy.linspace(-m, m, bin_num + 1)
        bins[-1] = bins[-1] * bound_scale
        bins[0] = bins[0] * bound_scale
    return bins


def get_chisq_grid(hist2d, grid_x, grid_y, G_PDF_bin, gh, bin_num, bin_num2):
    nums = numpy.zeros((bin_num,))

    for i in range(bin_num):
        x1 = 1. / gh * (grid_x - G_PDF_bin[i])
        x2 = 1. / gh * (grid_x - G_PDF_bin[i + 1])
        if gh > 0:
            idx1 = x1 > grid_y
            idx2 = x2 <= grid_y
        if gh < 0:
            idx1 = x1 <= grid_y
            idx2 = x2 > grid_y
        idx = idx1 & idx2
        nums[i] = numpy.sum(hist2d[idx])
    n1 = numpy.flip(nums[0:bin_num2],axis=0)
    n2 = nums[bin_num2:]

    return numpy.sum((n1 - n2) ** 2 / (n1 + n2)) * 0.5, n1, n2


def find_shear_grid(G, NU, G_PDF_bin, G_hist_bin, NU_hist_bin, left=-0.11, right=0.1,
                    chisq_gap=100, dg=0.0005, max_iters=40, fit_num=10, ax=False):
    data_num = G.shape[0]
    G = numpy.ascontiguousarray(G, dtype=numpy.float64)
    NU = numpy.ascontiguousarray(NU, dtype=numpy.float64)

    bin_num = G_PDF_bin.shape[0] - 1
    bin_num2 = int(bin_num * 0.5)
    inverse = range(int(bin_num / 2 - 1), -1, -1)

    num_g = G_hist_bin.shape[0] - 1
    num_nu = NU_hist_bin.shape[0] - 1
    grid_x, grid_y = numpy.zeros((num_nu, num_g), dtype=numpy.float64), numpy.zeros((num_nu, num_g),
                                                                                    dtype=numpy.float64)
    hist_num2d = numpy.zeros((num_nu, num_g), dtype=numpy.intc)

    hist2d_fast(G, NU, data_num, G_hist_bin, NU_hist_bin, num_g, num_nu, hist_num2d, grid_x, grid_y)
    idx = hist_num2d > 0
    hist_num2d = hist_num2d[idx]
    grid_x = grid_x[idx]
    grid_y = grid_y[idx]
    iters = 0
    change = 1

    while change == 1:
        change = 0
        mc = (left + right) / 2.
        mcl = left
        mcr = right
        fmc = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mc, bin_num, bin_num2)[0]
        fmcl = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mcl, bin_num, bin_num2)[0]
        fmcr = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mcr, bin_num, bin_num2)[0]
        #         print("%d. %.4f %.2f  %.4f %.2f  %.4f %.2f"%(iters, left, fmcl, mc,fmc,right,fmcr))
        temp = fmc + chisq_gap

        if fmcl > temp:
            #             left = (mc + mcl) / 2.
            left = mcl + (mc - mcl) / 3
            change = 1
        if fmcr > temp:
            #             right = (mc + mcr) / 2.
            right = mcr - (mcr - mc) / 3
            change = 1

        iters += 1
        if right - left < dg:
            break
        if iters > max_iters:
            break

    ghs = numpy.linspace(left, right, fit_num)
    xi2 = numpy.array([get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, gh, bin_num, bin_num2)[0] for gh in ghs])

    coeff = tool_box.fit_1d(ghs, xi2, 2, "scipy")
    gh = -coeff[1] / 2. / coeff[2]
    gh_sig = 0.70710678118 / numpy.sqrt(coeff[2])

    n1, n2 = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, gh, bin_num, bin_num2)[1:3]
    asym = (n1 ** 2 + n2 ** 2) / 2 / n1 / n2
    asym = numpy.sum(asym) - bin_num2

    chisqs_min = coeff[0] - coeff[1] ** 2 / 4 / coeff[2]
    chisqs_min_ = numpy.sum((n1-n2)**2/(n1+n2))*0.5
    if ax:
        ax.scatter(ghs, xi2)
        ax.plot(ghs, coeff[0] + coeff[1] * ghs + coeff[2] * ghs ** 2, c="C1")
        x1, x2 = left - (right - left)*0.2, right + (right - left)*0.2
        ax.plot([x1, x2], [chisqs_min, chisqs_min], ls="--", c="k", label="%.2f(%.2f)" % (chisqs_min,chisqs_min_))
        ax.set_xlim((x1,x2))
        ax.legend(loc="lower left")
        ax.set_title("asym=%.3e" % asym)
    return gh, gh_sig, coeff, asym, chisqs_min_, chisqs_min, ghs, xi2, grid_x, grid_y, hist_num2d


def get_chisq_grid_corr_new(hist2d, hist2d_corr, grid_x, grid_y, grid_x_corr, grid_y_corr, G_PDF_bin, gh, bin_num,
                            bin_num2):
    nums = numpy.zeros((bin_num,))
    for i in range(bin_num):
        x1 = 1. / gh * (grid_x - G_PDF_bin[i])
        x2 = 1. / gh * (grid_x - G_PDF_bin[i + 1])
        if gh > 0:
            idx1 = x1 > grid_y
            idx2 = x2 <= grid_y
        if gh < 0:
            idx1 = x1 <= grid_y
            idx2 = x2 > grid_y
        idx = idx1 & idx2

        #         correction
        x1_corr = 1. / gh * (grid_x_corr - G_PDF_bin[i])
        x2_corr = 1. / gh * (grid_x_corr - G_PDF_bin[i + 1])
        if gh > 0:
            idx1 = x1_corr > grid_y_corr
            idx2 = x2_corr <= grid_y_corr
        if gh < 0:
            idx1 = x1_corr <= grid_y_corr
            idx2 = x2_corr > grid_y_corr
        idx_corr = idx1 & idx2

        nums[i] = numpy.sum(hist2d[idx]) - numpy.sum(hist2d_corr[idx_corr])

    n1 = numpy.flip(nums[0:bin_num2],axis=0)
    n2 = nums[bin_num2:]

    return numpy.sum((n1 - n2) ** 2 / (n1 + n2)) * 0.5, n1, n2


def find_shear_grid_corr_new(G, NU, G_corr, NU_corr, G_PDF_bin, G_hist_bin, NU_hist_bin,
                             left=-0.11, right=0.1,chisq_gap=100, dg=0.005, max_iters=40, fit_num=10, ax=False):
    data_num = G.shape[0]
    data_num_corr = G_corr.shape[0]

    G = numpy.ascontiguousarray(G, dtype=numpy.float64)
    NU = numpy.ascontiguousarray(NU, dtype=numpy.float64)
    G_corr = numpy.ascontiguousarray(G_corr, dtype=numpy.float64)
    NU_corr = numpy.ascontiguousarray(NU_corr, dtype=numpy.float64)

    bin_num = G_PDF_bin.shape[0] - 1
    bin_num2 = int(bin_num * 0.5)
    inverse = range(int(bin_num / 2 - 1), -1, -1)

    num_g = G_hist_bin.shape[0] - 1
    num_nu = NU_hist_bin.shape[0] - 1
    grid_x, grid_y = numpy.zeros((num_nu, num_g), dtype=numpy.float64), numpy.zeros((num_nu, num_g),
                                                                                    dtype=numpy.float64)
    hist_num2d = numpy.zeros((num_nu, num_g), dtype=numpy.intc)

    grid_x_corr, grid_y_corr = numpy.zeros((num_nu, num_g), dtype=numpy.float64), numpy.zeros((num_nu, num_g),
                                                                                              dtype=numpy.float64)
    hist_num2d_corr = numpy.zeros((num_nu, num_g), dtype=numpy.intc)

    hist2d_fast(G, NU, data_num, G_hist_bin, NU_hist_bin, num_g, num_nu, hist_num2d, grid_x, grid_y)
    hist2d_fast(G_corr, NU_corr, data_num_corr, G_hist_bin, NU_hist_bin, num_g, num_nu, hist_num2d_corr, grid_x_corr,
                grid_y_corr)

    hist_all = hist_num2d + hist_num2d_corr
    idx = hist_all > 0
    hist_num2d = hist_num2d[idx]
    grid_x = grid_x[idx]
    grid_y = grid_y[idx]
    hist_num2d_corr = hist_num2d_corr[idx]
    grid_x_corr = grid_x_corr[idx]
    grid_y_corr = grid_y_corr[idx]

    iters = 0
    change = 1
    while change == 1:
        change = 0
        mc = (left + right) / 2.
        mcl = left
        mcr = right
        fmc = get_chisq_grid_corr_new(hist_num2d, hist_num2d_corr, grid_x, grid_y, grid_x_corr, grid_y_corr, G_PDF_bin,
                                      mc, bin_num, bin_num2)[0]
        fmcl = get_chisq_grid_corr_new(hist_num2d, hist_num2d_corr, grid_x, grid_y, grid_x_corr, grid_y_corr, G_PDF_bin,
                                       mcl, bin_num, bin_num2)[0]
        fmcr = get_chisq_grid_corr_new(hist_num2d, hist_num2d_corr, grid_x, grid_y, grid_x_corr, grid_y_corr, G_PDF_bin,
                                       mcr, bin_num, bin_num2)[0]
        #         print("%d. %.4f %.2f  %.4f %.2f  %.4f %.2f"%(iters, left, fmcl, mc,fmc,right,fmcr))
        temp = fmc + chisq_gap

        if fmcl > temp:
            #             left = (mc + mcl) / 2.
            left = mcl + (mc - mcl) / 3
            change = 1
        if fmcr > temp:
            #             right = (mc + mcr) / 2.
            right = mcr - (mcr - mc) / 3
            change = 1

        iters += 1
        if right - left < dg:
            break
        if iters > max_iters:
            break

    ghs = numpy.linspace(left, right, fit_num)
    xi2 = numpy.array([get_chisq_grid_corr_new(hist_num2d, hist_num2d_corr, grid_x, grid_y, grid_x_corr, grid_y_corr,
                                               G_PDF_bin, gh, bin_num, bin_num2)[:1] for gh in ghs])

    coeff = tool_box.fit_1d(ghs, xi2, 2, "scipy")
    chisqs_min = coeff[0] - coeff[1] ** 2 / 4 / coeff[2]
    gh = -coeff[1] / 2. / coeff[2]
    gh_sig = 0.70710678118 / numpy.sqrt(coeff[2])

    n1, n2 = get_chisq_grid_corr_new(hist_num2d, hist_num2d_corr, grid_x, grid_y, grid_x_corr, grid_y_corr, G_PDF_bin,
                                     gh, bin_num, bin_num2)[1:3]

    chisqs_min_ = numpy.sum((n1-n2)**2/(n1 + n2))*0.5
    asym = (n1 ** 2 + n2 ** 2) / 2 / n1 / n2
    asym = numpy.sum(asym) - bin_num2


    if ax:
        ax.scatter(ghs, xi2)
        ax.plot(ghs, coeff[0] + coeff[1] * ghs + coeff[2] * ghs ** 2, c="C1")
        x1, x2 = left - (right - left)*0.2, right + (right - left)*0.2
        ax.plot([x1, x2], [chisqs_min, chisqs_min], ls="--", c="k", label="%.2f(%.2f)" % (chisqs_min, chisqs_min_))
        ax.set_xlim((x1,x2))
        ax.legend(loc="lower left")
        ax.set_title("asym=%.3e"%asym)

    return gh, gh_sig, coeff, asym, chisqs_min_, chisqs_min, ghs, xi2, grid_x, grid_y, hist_num2d