from sys import path, argv
path.append("/home/hklee/work/mylib")
import time
from plot_tool import Image_Plot
import tool_box
from Fourier_Quad import Fourier_Quad
import numpy
import h5py
import time
from mpi4py import MPI



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


def find_shear_grid(G, NU, G_PDF_bin, G_hist_bin, NU_hist_bin, chisq_gap=100, dg=0.01):
    bin_num = G_PDF_bin.shape[0] - 1
    bin_num2 = int(bin_num * 0.5)
    inverse = range(int(bin_num / 2 - 1), -1, -1)

    num_g = G_hist_bin.shape[0] - 1
    num_nu = NU_hist_bin.shape[0] - 1
    grid_x, grid_y = numpy.zeros((num_nu, num_g)), numpy.zeros((num_nu, num_g))

    for i in range(num_g):
        grid_x[:, i] = (G_hist_bin[i] + G_hist_bin[i + 1]) / 2
    for i in range(num_nu):
        grid_y[i] = (NU_hist_bin[i] + NU_hist_bin[i + 1]) / 2

    hist2d = numpy.histogram2d(NU, G, [NU_hist_bin, G_hist_bin])[0]

    gl, gr = -0.1, 0.1

    num_gh = 100
    xi2 = numpy.zeros((num_gh,))
    ghs = numpy.linspace(gl, gr, num_gh)

    nums = numpy.zeros((bin_num,))

    #     for ig in range(num_gh):
    #         gh = ghs[ig]
    #         xi2[ig] = get_chisq_grid(hist2d, grid_x, grid_y,G_PDF_bin,gh, bin_num, bin_num2,inverse)

    fit_num = 20
    iters = 0
    change = 1
    left, right = -0.1, 0.099
    while change == 1:
        change = 0
        mc = (left + right) / 2.
        mcl = left
        mcr = right
        fmc = get_chisq_grid(hist2d, grid_x, grid_y, G_PDF_bin, mc, bin_num, bin_num2, inverse)
        fmcl = get_chisq_grid(hist2d, grid_x, grid_y, G_PDF_bin, mcl, bin_num, bin_num2, inverse)
        fmcr = get_chisq_grid(hist2d, grid_x, grid_y, G_PDF_bin, mcr, bin_num, bin_num2, inverse)

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
    xi2 = numpy.array([get_chisq_grid(hist2d, grid_x, grid_y, G_PDF_bin, gh, bin_num, bin_num2, inverse) for gh in ghs])
    print(iters, ghs)

    # img = Image_Plot()
    # img.subplots(1, 1)
    # img.axs[0][0].plot(ghs, xi2)
    # img.show_img()

    coeff = tool_box.fit_1d(ghs, xi2, 2, "scipy")
    gh = -coeff[1] / 2. / coeff[2]
    gh_sig = 0.70710678118 / numpy.sqrt(coeff[2])

    return gh, gh_sig, grid_x, grid_y, hist2d


data_path = argv[1]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

num_g, num_nu = 1500, 200

fq = Fourier_Quad(12,124)

data_type = ["noise_free", "noisy_cpp"]
result = []
for i in range(2):
    h5f = h5py.File(data_path + "/data_%d_%s.hdf5"%(rank, data_type[i]),"r")
    data = h5f["/data"][()]
    h5f.close()

    G1,G2 = data[:,0], data[:,1]
    NU1, NU2 = data[:,2]+data[:,3], data[:,2]-data[:,3]

    t1 = time.time()
    gh1, gh1_sig = fq.find_shear(G1, NU1, 10)[:2]

    t2 = time.time()
    gh2, gh2_sig = fq.find_shear(G2, NU2, 10)[:2]

    result.extend([gh1, gh1_sig, gh2, gh2_sig])


    t3 = time.time()

    G1_bin = tool_box.set_bin(G1, 10,100)
    G1_hist_bin = tool_box.set_bin(G1, num_g+1, 1.1)
    NU1_hist_bin = set_bin(NU1, num_nu+1, 1.1,"linear")

    gh1_grid, gh1_sig_grid, hist1, grid_G1, grid_NU1 = find_shear_grid(G1, NU1, G1_bin, G1_hist_bin, NU1_hist_bin,
                                                             chisq_gap=200,dg=0.005)


    t4 = time.time()

    G2_bin = tool_box.set_bin(G2, 10,100)
    G2_hist_bin = tool_box.set_bin(G2, num_g+1, 1.1)
    NU2_hist_bin = set_bin(NU2, num_nu+1, 1.1,"linear")

    gh2_grid, gh2_sig_grid, hist2, grid_G2, grid_NU2 = find_shear_grid(G2, NU2, G2_bin, G2_hist_bin, NU2_hist_bin,
                                                             chisq_gap=500,dg=0.005)
    result.extend([gh1_grid, gh1_sig_grid, gh2_grid, gh2_sig_grid])
    t5 = time.time()

    print("RANK %d: %.2f, %.2f, %.2f, %.2f\n"%(rank,t2-t1, t3-t2, t4-t3, t5-t4))

comm.Barrier()

results = comm.gather(result, root=0)

if rank == 0:
    numpy.savez(data_path + "/cache.npz", numpy.array(results))

comm.Barrier()