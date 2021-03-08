from sys import path,argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
from Fourier_Quad import Fourier_Quad
import numpy
import h5py
import time
import ctypes
import numpy.ctypeslib as ctl
from mpi4py import MPI

histlib = ctypes.cdll.LoadLibrary("./libhist.so")

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


def set_bin(data, bin_num, bound_scale, method="pmean"):
    if method == "pmean":
        temp_data = numpy.sort(numpy.abs(data))
        bin_size = len(temp_data) / bin_num
        bins = numpy.array([temp_data[int(i * bin_size)] for i in range(bin_num)])
        bins[-1] = bins[-1] * bound_scale
    elif method == "linear":
        m1, m2 = data.min(), data.max()
        bins = numpy.linspace(m1, m2, bin_num + 1)
        bins[-1] = bins[-1] * bound_scale
    else:
        m = max(numpy.abs(data.min()), numpy.abs(data.max()))
        bins = numpy.linspace(-m, m, bin_num + 1)
        bins[-1] = bins[-1] * bound_scale
        bins[0] = bins[0] * bound_scale
    return bins

def set_bin_(data, bin_num, bound_scale, method="mean"):
    if method == "mean":
        temp_data = numpy.sort(data)
        bin_size = len(temp_data) / bin_num
        bins = [temp_data[int(i * bin_size)] for i in range(bin_num)]
        bins.append(temp_data[-1] * bound_scale)
        bins = numpy.array(bins)

    elif method == "linear":
        m1, m2 = data.min(), data.max()
        bins = numpy.linspace(m1, m2, bin_num + 1)
        bins[-1] = bins[-1] * bound_scale
    elif method == "log":
        data_min = data.min()
        bin_num2 = int(bin_num / 2)
        bins = numpy.zeros((bin_num + 1,))

        if data_min < 0:
            temp_data = numpy.sort(numpy.abs(data))
            data_min, data_max = temp_data[0] * 0.98, temp_data[-1]

            bin_num_ = bin_num2 - 1
            invers = range(bin_num_, -1, -1)
            hbins = tool_box.set_bin_log(data_min, data_max, bin_num_)

            bins[bin_num2 + 1:] = hbins
            bins[:bin_num2] = -hbins[invers]
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


def find_shear_grid(G, NU, G_PDF_bin, G_hist_bin, NU_hist_bin, chisq_gap=100, dg=0.01, fit_num=10):
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
    left, right = -0.1, 0.099
    while change == 1:
        change = 0
        mc = (left + right) / 2.
        mcl = left
        mcr = right
        fmc = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mc, bin_num, bin_num2, inverse)
        fmcl = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mcl, bin_num, bin_num2, inverse)
        fmcr = get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, mcr, bin_num, bin_num2, inverse)

        temp = fmc + chisq_gap

        if fmcl > temp:
            left = (mc + mcl) / 2.
            change = 1
        if fmcr > temp:
            right = (mc + mcr) / 2.
            change = 1
        #         print(iters, left,right)
        iters += 1
        if right - left < dg:
            break
        if iters > 12:
            break

    ghs = numpy.linspace(left, right, fit_num)
    xi2 = numpy.array(
        [get_chisq_grid(hist_num2d, grid_x, grid_y, G_PDF_bin, gh, bin_num, bin_num2, inverse) for gh in ghs])
    #     print(iters, ghs)

    #     img = Image_Plot()
    #     img.subplots(1, 1)
    #     img.axs[0][0].plot(ghs, xi2)
    #     img.show_img()

    coeff = tool_box.fit_1d(ghs, xi2, 2, "scipy")
    gh = -coeff[1] / 2. / coeff[2]
    gh_sig = 0.70710678118 / numpy.sqrt(coeff[2])

    return gh, gh_sig, grid_x, grid_y, hist_num2d



data_path = argv[1]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

num_g, num_nu = 500, 500

fq = Fourier_Quad(12,124)

data_type = argv[2]#["noise_free", "noisy_cpp"]
result = []


h5f = h5py.File(data_path + "/data_%d_%s.hdf5"%(rank, data_type),"r")
data = h5f["/data"][()]
h5f.close()

G1 = numpy.ascontiguousarray(data[:, 0], dtype=numpy.float64)
G2 = numpy.ascontiguousarray(data[:, 1], dtype=numpy.float64)
NU1 = numpy.ascontiguousarray(data[:, 2] + data[:, 3], dtype=numpy.float64)
NU2 = numpy.ascontiguousarray(data[:, 2] - data[:, 3], dtype=numpy.float64)

G1_bin = tool_box.set_bin(G1, 10, 100)
# G1_hist_bin = tool_box.set_bin(G1, num_g, 1.001)
G1_hist_bin = set_bin_(G1, num_g, 1.001,"log")
NU1_hist_bin = set_bin_(NU1, num_nu, 1.001, "log")

G2_bin = tool_box.set_bin(G2, 10, 100)
# G2_hist_bin = tool_box.set_bin(G2, num_g, 1.001,"log")
G2_hist_bin = set_bin_(G2, num_g, 1.001,"log")
NU2_hist_bin = set_bin_(NU2, num_nu, 1.001, "log")

print(NU1.min(), NU2.max(), NU2.min(), NU2.max())
t1 = time.time()

gh1, gh1_sig = fq.find_shear(G1, NU1, 10)[:2]

t2 = time.time()
gh2, gh2_sig = fq.find_shear(G2, NU2, 10)[:2]

result.extend([gh1, gh1_sig, gh2, gh2_sig])

t3 = time.time()

gh1_grid, gh1_sig_grid, hist1, grid_G1, grid_NU1 = find_shear_grid(G1, NU1, G1_bin, G1_hist_bin, NU1_hist_bin,
                                                         chisq_gap=200,dg=0.005)

t4 = time.time()

gh2_grid, gh2_sig_grid, hist2, grid_G2, grid_NU2 = find_shear_grid(G2, NU2, G2_bin, G2_hist_bin, NU2_hist_bin,
                                                         chisq_gap=200,dg=0.005)
result.extend([gh1_grid, gh1_sig_grid, gh2_grid, gh2_sig_grid])
t5 = time.time()

print("RANK %d: %.2f, %.2f, %.2f, %.2f\n"%(rank,t2-t1, t3-t2, t4-t3, t5-t4))

comm.Barrier()

results = comm.gather(result, root=0)


if rank == 0:

    result_arr = numpy.array(results)

    h5f = h5py.File(data_path + "/shear.hdf5","r")
    g1 = h5f["/g1"][()]
    g2 = h5f["/g2"][()]
    h5f.close()

    mc1 = tool_box.data_fit(g1, result_arr[:,0], result_arr[:,1])
    mc2 = tool_box.data_fit(g2, result_arr[:,2], result_arr[:,3])

    mc1_grid = tool_box.data_fit(g1, result_arr[:,4], result_arr[:,5])
    mc2_grid = tool_box.data_fit(g2, result_arr[:,6], result_arr[:,7])

    numpy.savez(data_path + "/cache_%s.npz"%data_type, result_arr, g1, mc1, mc1_grid, g2, mc2, mc2_grid)

comm.Barrier()