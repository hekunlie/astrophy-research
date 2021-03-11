from sys import path, argv
path.append("/home/hklee/work/mylib")
import time
from plot_tool import Image_Plot
import tool_box
from Fourier_Quad import Fourier_Quad
import numpy
import h5py
import time
import ctypes
import numpy.ctypeslib as ctl
from mpi4py import MPI


histlib = ctypes.cdll.LoadLibrary("/home/hklee/work/mylib/c4py.so")

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


def set_bin(data, bin_num, bound_scale, method="mean", log_end = 5):
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
    elif method == "linear":
        m1, m2 = data.min(), data.max()
        bins = numpy.linspace(m1, m2, bin_num + 1)
        bins[-1] = bins[-1] * bound_scale

    elif method == "mix_bin":
        mean_num, log_num = bin_num
        mean_num2, log_num2 = int(mean_num/2), int(log_num/2)

        hbins = numpy.zeros((mean_num2 + log_num2,))

        temp_data = numpy.sort(numpy.abs(data))  # [:int(len(data[data>0])*0.99)]

        bin_size = len(temp_data) / mean_num * 2
        mean_bins = numpy.array([temp_data[int(i * bin_size)] for i in range(1, mean_num2)])

        hbins[:mean_num2-1] = mean_bins

        hbins[mean_num2-1] = temp_data[-1]

        ed = numpy.log10(temp_data[-1])
        log_bin = 10**numpy.linspace(ed, ed + log_end, log_num2+1)
        hbins[mean_num2:] = log_bin[1:]

        inv = range(mean_num2+ log_num2-1, -1, -1)
        bins = numpy.zeros((mean_num+ log_num+1))
        bins[mean_num2+ log_num2+1:] = hbins
        bins[:mean_num2+log_num2] = -hbins[inv]

    elif method == "log":
        data_min = data.min()
        bin_num2 = int(bin_num / 2)
        bins = numpy.zeros((bin_num + 1,))

        if data_min < 0:
            temp_data = numpy.sort(numpy.abs(data))
            data_min, data_max = temp_data[0], temp_data[-1]
            if data_min < 10:
                data_min = 10
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
            bins = tool_box.set_bin_log(data_min, data_max, bin_num)
            bins[0] = bins[0] * 0.95
            bins[-1] = bins[-1] * bound_scale
    else:
        m = max(numpy.abs(data.min()), numpy.abs(data.max()))
        bins = numpy.linspace(-m, m, bin_num + 1)
        bins[-1] = bins[-1] * bound_scale
        bins[0] = bins[0] * bound_scale
    return bins


def get_chisq_grid(hist2d, grid_x, grid_y, G_PDF_bin, gh, bin_num, bin_num2, inverse):
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
    n1 = nums[0:bin_num2][inverse]
    n2 = nums[bin_num2:]

    return numpy.sum((n1 - n2) ** 2 / (n1 + n2)) * 0.5


def find_shear_grid(G, NU, G_PDF_bin, G_hist_bin, NU_hist_bin, chisq_gap=100, dg=0.01, fit_num=10, ax=False):
    data_num = G.shape[0]

    bin_num = G_PDF_bin.shape[0] - 1
    bin_num2 = int(bin_num * 0.5)
    inverse = range(int(bin_num / 2 - 1), -1, -1)

    num_g = G_hist_bin.shape[0] - 1
    num_nu = NU_hist_bin.shape[0] - 1
    grid_x, grid_y = numpy.zeros((num_nu, num_g), dtype=numpy.float64), numpy.zeros((num_nu, num_g),
                                                                                    dtype=numpy.float64)
    hist_num2d = numpy.zeros((num_nu, num_g), dtype=numpy.intc)

    hist2d_fast(G, NU, data_num, G_hist_bin, NU_hist_bin, num_g, num_nu, hist_num2d, grid_x, grid_y)

    iters = 0
    change = 1
    left, right = -0.1, 0.11
    while change == 1:
        change = 0
        mc = (left + right) / 2.
        mcl = left
        mcr = right
        fmc = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mc, bin_num, bin_num2, inverse)
        fmcl = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mcl, bin_num, bin_num2, inverse)
        fmcr = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mcr, bin_num, bin_num2, inverse)
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
        if iters > 40:
            break

    ghs = numpy.linspace(left, right, fit_num)
    xi2 = numpy.array(
        [get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, gh, bin_num, bin_num2, inverse) for gh in ghs])
    #     print(iters, ghs)
    if ax:
        ax.plot(ghs, xi2)

    coeff = tool_box.fit_1d(ghs, xi2, 2, "scipy")
    gh = -coeff[1] / 2. / coeff[2]
    gh_sig = 0.70710678118 / numpy.sqrt(coeff[2])

    return gh, gh_sig, grid_x, grid_y, hist_num2d


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

itemsize = MPI.DOUBLE.Get_size()
element_num = 30*5*4
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
result = numpy.ndarray(buffer=buf1, dtype='d', shape=(5,4,30)) # array filled with zero


bins = argv[1]

data_path = "/mnt/perc/hklee/PDF_test/PDF_grid_test/data"
h5f = h5py.File(data_path +"/shear.hdf5","r")
g1 = h5f["/g1"][()]
g2 = h5f["/g2"][()]
h5f.close()

shear_tag = rank
print(shear_tag, g1[shear_tag], g2[shear_tag])

h5f = h5py.File(data_path + "/data_%d_noisy_cpp.hdf5" % shear_tag, "r")
data = h5f["/data"][()]
h5f.close()

G1 = numpy.ascontiguousarray(data[:, 0], dtype=numpy.float64)
G2 = numpy.ascontiguousarray(data[:, 1], dtype=numpy.float64)
NU1 = numpy.ascontiguousarray(data[:, 2] + data[:, 3], dtype=numpy.float64)
NU2 = numpy.ascontiguousarray(data[:, 2] - data[:, 3], dtype=numpy.float64)

G1_bin = tool_box.set_bin(G1, 10, 100)
G2_bin = tool_box.set_bin(G2, 10, 100)

img = Image_Plot(xpad=0.2, ypad=0.2)
img.subplots(2, 5)

bin_type = [bins for i in range(5)]

bins_nums = [[100, 50], [100, 100], [200, 100], [200, 200], [500, 200]]

for j in range(5):
    num_g, num_nu = bins_nums[j]

    G1_hist_bin = set_bin(G1[:3000], bins_nums[j], 1.001, bin_type[j])
    NU1_hist_bin = set_bin(NU1[:3000], bins_nums[j], 1.001, bin_type[j])

    G2_hist_bin = set_bin(G2[:3000], bins_nums[j], 1.001, bin_type[j])
    NU2_hist_bin = set_bin(NU2[:3000], bins_nums[j], 1.001, bin_type[j])

    # G1_hist_bin = set_bin(G1[:3000], num_g, 1.001, bin_type[j])
    # NU1_hist_bin = set_bin(NU1[:3000], num_nu, 1.001, bin_type[j])
    #
    # G2_hist_bin = set_bin(G2[:3000], num_g, 1.001, bin_type[j])
    # NU2_hist_bin = set_bin(NU2[:3000], num_nu, 1.001, bin_type[j])

    t1 = time.time()
    gh1_grid, gh1_sig_grid, hist1, grid_G1, grid_NU1 = find_shear_grid(G1, NU1, G1_bin, G1_hist_bin, NU1_hist_bin,
                                                                       chisq_gap=100, dg=0.005, ax=img.axs[0][j])
    t2 = time.time()
    result_str = "g1: %.5f. measured: %.5f(%.5f).  diff: %.5f(%.5f). time: %.2f sec" % (
    g1[shear_tag], gh1_grid, gh1_sig_grid, g1[shear_tag] - gh1_grid, gh1_sig_grid, t2 - t1)
    print(result_str)

    t3 = time.time()
    gh2_grid, gh2_sig_grid, hist2, grid_G2, grid_NU2 = find_shear_grid(G2, NU2, G2_bin, G2_hist_bin, NU2_hist_bin,
                                                                       chisq_gap=100, dg=0.005, ax=img.axs[1][j])
    t4 = time.time()
    result_str = "g2: %.5f. measured: %.5f(%.5f).  diff: %.5f(%.5f). time: %.2f sec" % (
    g2[shear_tag], gh2_grid, gh2_sig_grid, g2[shear_tag] - gh2_grid, gh2_sig_grid, t4 - t3)
    print(result_str)

    result[j, :, rank] = gh1_grid, gh1_sig_grid, gh2_grid, gh2_sig_grid

    text = "%s [%d, %d]\ng1: %.5f,  %.5f(%.5f)\ndiff: %.5f" % (
    bin_type[j], num_g, num_nu, g1[shear_tag], gh1_grid, gh1_sig_grid, g1[shear_tag] - gh1_grid)
    img.axs_text(0, j, 0.9, 0.05, text, text_fontsize=15)
    text = "%s [%d, %d]\ng2: %.5f,  %.5f(%.5f)\ndiff: %.5f" % (
    bin_type[j], num_g, num_nu, g2[shear_tag], gh2_grid, gh2_sig_grid, g2[shear_tag] - gh2_grid)
    img.axs_text(1, j, 0.9, 0.05, text, text_fontsize=15)

img.save_img("./pic/%s_%d.png"%(bins,rank))
# img.show_img()
img.close_img()
comm.Barrier()

if rank == 0:
    img = Image_Plot()
    img.subplots(2, 5)
    for i in range(5):
        img.axs[0][i].errorbar(g1, result[i, 0,:], result[i, 1,:], marker="s", ms=5, fmt=" ", label="g1")
        img.axs[0][i].errorbar(g2, result[i, 2,:], result[i, 3,:], marker="o", ms=5, fmt=" ", label="g2")

        img.axs[1][i].errorbar(g1, (g1 - result[i, 0,:]) * 100, result[i, 1,:], marker="s", ms=5, fmt=" ", label="g1")
        img.axs[1][i].errorbar(g2, (g2 - result[i, 2,:]) * 100, result[i, 3,:], marker="o", ms=5, fmt=" ", label="g2")

        mc1 = tool_box.data_fit(g1, result[i, 0,:], result[i, 1,:])
        mc2 = tool_box.data_fit(g2, result[i, 2,:], result[i, 3,:])

        strs = "m1: %.4f(%.4f), c1: %.5f(%.5f)\nm2: %.4f(%.4f), c2: %.5f(%.5f)" % (
        mc1[0] - 1, mc1[1], mc1[2], mc1[3], mc2[0] - 1, mc2[1], mc2[2], mc2[3])
        img.axs_text(0, i, 0.95, 0.05, strs, text_fontsize=12)

        img.axs[0][i].plot([-0.04, 0.04], [-0.04, 0.04], ls="--", alpha=0.3, c="gray")
        img.axs[1][i].plot([-0.04, 0.04], [0, 0], ls="--", alpha=0.3, c="gray")

        img.axs[0][i].legend(loc="lower right")
        img.axs[1][i].legend(loc="lower right")
        img.set_label(0, i, 0, "Measured g")
        img.set_label(0, i, 1, "True g")

        img.set_label(1, i, 0, "diff *100")
        img.set_label(1, i, 1, "True g")

        print(mc1)
        print(mc2)
    img.save_img("./pic/%s_mc.png" % bins)
    # img.show_img()
    img.close_img()

comm.Barrier()