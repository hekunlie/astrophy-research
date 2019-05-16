import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import time
from Fourier_Quad import Fourier_Quad
import tool_box
import plot_tool
import h5py
from mpi4py import MPI
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

area_id = int(argv[1])
foreground = argv[2]

itemsize = MPI.DOUBLE.Get_size()
element_num = 13*4
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf1, itemsize = win1.Shared_query(0)


# create a numpy array from buf
# code can run successfully without the following step
# buf = np.array(buf, dtype='float64', copy=False) # may be redundant
# "d" means double = 'float64'
result = numpy.ndarray(buffer=buf1, dtype='d', shape=(4,13)) # array filled with zero


tag = rank


# data_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/w_%d/"%area_id
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/%s/w_%d/"%(foreground,area_id)
h5f = h5py.File(data_path + "%d.hdf5"%tag,"a")

print(list(h5f.keys()))
data = h5f["/data_0"].value
print(data.shape)

try:
    del h5f["/result"]
except:
    pass
gt = data[:, 0]
gx = data[:, 1]
n = data[:, 2]
u = data[:, 3]
crit = data[:, 4]*3.882833518*100
g1 = data[:, 5]

gt_c = gt*crit
gx_c = gx*crit
nu1 = n + u
nu2 = n - u

print(gt)
print(g1)
print(crit)

gtc_bin = tool_box.set_bin(gt_c,8, 100)
gt_bin = tool_box.set_bin(gt,8, 100)

fq = Fourier_Quad(1, 1)


g_corr_t, corr_sig_t = fq.fmin_g_new(g=gt_c, nu=nu1, bin_num=10, scale=100,
                                     pic_path="%s/%d_t.png"%(data_path,tag), left=-500, right=500)

g_corr_x, corr_sig_x = fq.fmin_g_new(g=gx_c, nu=nu2, bin_num=10, scale=100,
                                     pic_path="%s/%d_x.png"%(data_path,tag), left=-500, right=500)
print(g_corr_t,corr_sig_t)
print(g_corr_x,corr_sig_x)

result[0, rank] = g_corr_t
result[1, rank] = corr_sig_t
result[2, rank] = g_corr_x
result[3, rank] = corr_sig_x

h5f["/result"] = [g_corr_t,corr_sig_t,g_corr_x,corr_sig_x]
h5f.close()

chi_t = []
chi_x = []

# img = plot_tool.Image_Plot()
# img.plot_img(1, 1)

# gh = numpy.linspace(-200, 200, 50)
#
# for i in gh:
#     chi = fq.G_bin(gt_c, nu1, i,gtc_bin, 0)
#     chi_t.append(chi)
#     chi = fq.G_bin(gx_c, nu2, i, gtc_bin, 0)
#     chi_x.append(chi)
#
# chi_t = numpy.array(chi_t)
# chi_x = numpy.array(chi_x)
# img.axs[0][0].scatter(gh, chi_t, label="T")
# img.axs[0][0].scatter(gh, chi_x, label="X")
#
# img.axs[0][0].legend()
# img.save_img()

# coeff_t = tool_box.fit_1d(gh, chi_t, 2, "scipy")
# corr_sig_t = numpy.sqrt(1 / 2. / coeff_t[2])
# g_corr_t = -coeff_t[1] / 2. / coeff_t[2]
#
# coeff_x = tool_box.fit_1d(gh, chi_x, 2, "scipy")
# corr_sig_x = numpy.sqrt(1 / 2. / coeff_x[2])
# g_corr_x = -coeff_x[1] / 2. / coeff_x[2]

comm.Barrier()
# plt.close()
if rank == 0:
    numpy.savez(data_path + "result_%d.npz"%area_id, result)
    img = plot_tool.Image_Plot()
    img.create_subfig(1, 2)
    gh = 10**numpy.linspace(numpy.log10(0.04), numpy.log10(15), 13)
    img.axs[0][0].errorbar(gh, result[0], result[1], label="T",capsize=img.cap_size)
    img.axs[0][0].errorbar(gh, result[2], result[3], label="X",capsize=img.cap_size)
    img.axs[0][0].set_xscale("log")
    img.axs[0][0].legend()

    img.axs[0][1].errorbar(gh, numpy.abs(result[0]), result[1], label="T",capsize=img.cap_size)
    img.axs[0][1].errorbar(gh, numpy.abs(result[2]), result[3], label="X",capsize=img.cap_size)
    img.axs[0][1].set_xscale("log")
    img.axs[0][1].set_yscale("log")
    img.axs[0][1].set_ylim(0.01, 200)
    img.axs[0][1].legend()

    img.save_img(data_path + "result_%d.pdf"%area_id)
    print(result)

