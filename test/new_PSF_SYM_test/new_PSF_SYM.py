import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import tool_box
import numpy
from mpi4py import MPI
import h5py


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


def set_bin(data, bin_num, scale=100., sort_method="abs", sym=True):
    """
    set up bins for 1-D data
    :param data:
    :param bin_num: total number of bins
    :param sort_method: "abs": set the scale of bins according to the absolute value of data
                        "posi": set the scale of bins according to the positive value
                        else: just set scale of bins from the small end to the large end
    :param sym: True: set up bins symmetrically, False: the bins for the positive values
    :return: bins, (N, ) numpy array
    """
    if sort_method == "abs":
        temp_data = numpy.sort(numpy.abs(data))
    elif sort_method == "posi":
        temp_data = numpy.sort(data[data > 0])
    else:
        temp_data = numpy.sort(data)

    if sym:
        bin_size = len(temp_data) / bin_num * 2
        bins = numpy.array([temp_data[int(i * bin_size)] for i in range(1, int(bin_num / 2))])
        bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
        bound = numpy.max(numpy.abs(data)) * scale
        bins = numpy.append(-bound, numpy.append(bins, bound))
    else:
        bin_size = len(temp_data) / bin_num
        bins = numpy.array([temp_data[int(i * bin_size)] for i in range(1, bin_num)])
        if sort_method == "abs" or sort_method == "posi":
            bins = numpy.sort(numpy.append([0], bins))
            bound = numpy.max(numpy.abs(data)) * scale
            bins = numpy.append(bins, bound)
        else:
            # for the sort of negative data
            bound = numpy.min(data) * scale
            bins = numpy.append(bound, numpy.append(bins, -bound))
    return bins

def chisq(g, nu, g_h, bins, bin_num2, inverse, ig_num):  # checked 2017-7-9!!!
    r"""
    to calculate the symmetry the shear estimators
    :param g: estimators from Fourier quad, 1-D numpy array
    :param nu: N + U for g1, N - U for g2, 1-D numpy array
    :param g_h: pseudo shear (guess)
    :param bins: bin of g for calculation of the symmetry, 1-D numpy array
    :param ig_num: the number of inner grid of bin to be neglected
    :return: chi square
    """
    G_h = g - nu * g_h
    num = numpy.histogram(G_h, bins)[0]
    n1 = num[0:bin_num2][inverse]
    n2 = num[bin_num2:]
    xi = (n1 - n2) ** 2 / (n1 + n2)
    return numpy.sum(xi[:len(xi) - ig_num]) * 0.5

def chisq_new(g, nu, g_h, bins, bin_num2, inverse, ig_num, num_exp):  # checked 2017-7-9!!!
    r"""
    to calculate the symmetry the shear estimators
    :param g: estimators from Fourier quad, 1-D numpy array
    :param nu: N + U for g1, N - U for g2, 1-D numpy array
    :param g_h: pseudo shear (guess)
    :param bins: bin of g for calculation of the symmetry, 1-D numpy array
    :param ig_num: the number of inner grid of bin to be neglected
    :return: chi square
    """
    G_h = g - nu * g_h
    num = numpy.histogram(G_h, bins)[0]
    n1 = num[0:bin_num2][inverse]
    n2 = num[bin_num2:]
    xi = ((n1 - num_exp)**2 + (n2 - num_exp)**2)/(n1 + n2)
    return numpy.sum(xi[:len(xi) - ig_num]) * 0.5

def shear_fit(x, y, fig_ax,label):
    coeff = tool_box.fit_1d(x, y, 2, "scipy")
    # y = a1 + a2*x a3*x^2 = a2(x+a1/2/a2)^2 +...
    # gh = - a1/2/a2, gh_sig = \sqrt(1/2/a2)
    g_h = -coeff[1] / 2. / coeff[2]
    g_sig = 0.70710678118 / numpy.sqrt(coeff[2])
    if fig_ax:
        fig_ax.scatter(x, y, alpha=0.7, s=5)
        fig_ax.plot(x, coeff[0] + coeff[1] * x + coeff[2] * x ** 2, alpha=0.7, label=label)
    return g_h, g_sig, coeff


def get_exp_num(mg, bin_num):
    bin_num2 = int(bin_num * 0.5)
    bins = set_bin(mg, bin_num, 10)

    inverse = range(int(bin_num2 - 1), -1, -1)
    num_ini = numpy.histogram(mg, bins)[0]
    n1 = num_ini[0:bin_num2][inverse]
    n2 = num_ini[bin_num2:]
    num_exp = (n1 + n2) / 2
    return num_exp, inverse, bins


result_row, result_col = 8, 10
itemsize = MPI.DOUBLE.Get_size()
element_num = result_row*result_col
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
result = numpy.ndarray(buffer=buf1, dtype='d', shape=(result_row,result_col)) # array filled with zero

fq = Fourier_Quad(12,13)

# # normal distribution
# num = 11
# shear = numpy.linspace(-0.05, 0.05, num)
# result = numpy.zeros((4, num))
# for i in range(num):
#
#     noise = 0 #numpy.random.normal(0, 0.01, 100000)
#     mg = numpy.random.normal(shear[i], 0.3, 10000) + noise
#     mn = numpy.ones_like(mg)
#
#     # img = Image_Plot()
#     # img.subplots(1, 2)
#     # img.axs[0][0].hist(mg, 100)
#     # gh, gh_sig = fmin_g_new(mg, mn, 8, fig_ax=img.axs[0][1])[:2]
#     gh1, gh1_sig = fq.find_shear(mg, mn, 8)[:2]
#     gh2, gh2_sig = fq.find_shear_new(mg, mn, 8)[:2]
#
#     result[0, i] = gh1
#     result[1, i] = gh1_sig
#     result[2, i] = gh2
#     result[3, i] = gh2_sig
#     print("Est g: %.5f (%.6f) %.5f (%.6f)"%(gh1, gh1_sig,gh2, gh2_sig))
#
#     # img.show_img()
#
# img = Image_Plot()
# img.subplots(1, 1)
# img.axs[0][0].errorbar(shear, result[0], result[1], capsize=4, marker="s", mfc="none",c="C2",label="Original $\chi^2$")
# img.axs[0][0].errorbar(shear, result[2], result[3], capsize=4, marker="s", mfc="none",c="C1",label="New $\chi^2$")
# img.axs[0][0].plot([-0.06, 0.06], [-0.06, 0.06], linestyle="--", c="grey")
# img.set_label(0,0,0,"$g_m$")
# img.set_label(0,0,1,"$g_t$")
# img.axs[0][0].legend(fontsize=img.legend_size)
# img.save_img("normal_pts.png")
# img.show_img()
#
# noise = numpy.random.normal(0, 0.03, 50000)
# mg = numpy.random.normal(shear[6], 0.3, 50000) + noise
# mn = numpy.ones_like(mg)
#
# gh1, gh1_sig = fq.find_shear(mg, mn, 8)[:2]
# gh2, gh2_sig = fq.find_shear_new(mg, mn, 8)[:2]
# print("Est g: %.5f (%.6f) %.5f (%.6f)"%(gh1, gh1_sig,gh2, gh2_sig))
#
# bin_num = 14
# bin_num2 = int(bin_num*0.5)
# scale = 10
# bins = set_bin(mg, bin_num, scale)
#
# inverse = range(int(bin_num2- 1), -1, -1)
# num_ini = numpy.histogram(mg, bins)[0]
# n1 = num_ini[0:bin_num2][inverse]
# n2 = num_ini[bin_num2:]
# num_exp = (n1 + n2) / 2
#
# fit_range = numpy.linspace(-0.1, 0.1, 200)
#
# chi_sq1 = numpy.array([chisq(mg, mn, g_hat, bins, bin_num2, inverse, 0) for g_hat in fit_range])
# chi_sq2 = numpy.array([chisq_new(mg, mn, g_hat, bins, bin_num2, inverse, 0, num_exp) for g_hat in fit_range])
#
# img = Image_Plot()
# img.subplots(1, 1)
# shear_fit(fit_range, chi_sq1, img.axs[0][0],"Original $\chi^2$")
# shear_fit(fit_range, chi_sq2, img.axs[0][0],"new $\chi^2$")
# img.set_label(0,0,0,"$\chi^2$")
# img.set_label(0,0,1,"$\hat g$")
# img.axs[0][0].legend(fontsize=img.legend_size)
# img.save_img("chi2.png")
# img.show_img()

data_path = "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/"
h5f = h5py.File(data_path + "result/data/data_%d.hdf5"%rank,"r")
data = h5f["/data"].value[:200000]
mg1 = data[:,2]
mg2 = data[:,3]
mnu1 = data[:,4] + data[:,5]
mnu2 = data[:,4] - data[:,5]

# gh1, gh1_sig = fq.find_shear(mg1, mnu1, 8)[:2]
# gh2, gh2_sig = fq.find_shear(mg2, mnu2, 8)[:2]
#
# gh1_new, gh1_sig_new = fq.find_shear_new(mg1, mnu1, 8)[:2]
# gh2_new, gh2_sig_new = fq.find_shear_new(mg2, mnu2, 8)[:2]
#
#
# result[0,rank] = gh1
# result[1,rank] = gh1_sig
# result[2,rank] = gh2
# result[3,rank] = gh2_sig
#
# result[4,rank] = gh1_new
# result[5,rank] = gh1_sig_new
# result[6,rank] = gh2_new
# result[7,rank] = gh2_sig_new


fit_range = numpy.linspace(-0.1, 0.1, 200)

bin_num = 8
bin_num2 = int(bin_num*0.5)

# G1
num_exp1, inverse1, bins1 = get_exp_num(mg1, bin_num)
# G2
num_exp2, inverse2, bins2 = get_exp_num(mg2, bin_num)

# original chi squared
chi_sq1 = numpy.array([chisq(mg1, mnu1, g_hat, bins1, bin_num2, inverse1, 0) for g_hat in fit_range])
chi_sq2 = numpy.array([chisq(mg2, mnu2, g_hat, bins2, bin_num2, inverse2, 0) for g_hat in fit_range])
# new chi squared
chi_sq1_new = numpy.array([chisq_new(mg1, mnu1, g_hat, bins1, bin_num2, inverse1, 0, num_exp1) for g_hat in fit_range])
chi_sq2_new = numpy.array([chisq_new(mg2, mnu2, g_hat, bins2, bin_num2, inverse2, 0, num_exp2) for g_hat in fit_range])

# plot
img = Image_Plot()
img.subplots(1, 2)

shear_fit(fit_range, chi_sq1, img.axs[0][0],"$\chi^2$ of $g_1$")
shear_fit(fit_range, chi_sq2, img.axs[0][0],"$\chi^2$ of $g_2$")

shear_fit(fit_range, chi_sq1_new, img.axs[0][1],"$\chi^2$ of $g_1$")
shear_fit(fit_range, chi_sq2_new, img.axs[0][1],"$\chi^2$ of $g_2$")

titles = ["Original $\chi^2$", "New $\chi^2$"]
for i in range(2):
    img.set_label(0,i,0,"$\chi^2$")
    img.set_label(0,i,1,"$\hat g$")
    img.axs[0][i].legend(fontsize=img.legend_size)
    img.axs[0][i].set_title(titles[i])

img.save_img("chisq/chisq_%d.png"%rank)


comm.Barrier()

# if rank == 0:
#
#     shear = numpy.loadtxt(data_path + "parameters/shear.dat")
#     g1 = shear[:result_col]
#     g2 = shear[result_col:]
#
#     img = Image_Plot()
#     img.subplots(1,2)
#
#     img.axs[0][0].errorbar(g1, result[0], result[1], fmt="none", capsize=4, marker="s", mfc="none",c="C2",label="Original $\chi^2$")
#     img.axs[0][0].scatter(g1, result[0], marker="s", c="C2")
#     img.axs[0][0].errorbar(g1, result[4], result[5], fmt="none", capsize=4, marker="s", mfc="none",c="C1",label="New $\chi^2$")
#     img.axs[0][0].scatter(g1, result[4], marker="s", c="C1")
#
#     img.axs[0][1].errorbar(g2, result[2], result[3], fmt="none",capsize=4, marker="s", mfc="none",c="C2",label="Original $\chi^2$")
#     img.axs[0][1].scatter(g2, result[2], marker="s", c="C2")
#     img.axs[0][1].errorbar(g2, result[6], result[7], fmt="none",capsize=4, marker="s", mfc="none",c="C1",label="New $\chi^2$")
#     img.axs[0][1].scatter(g2, result[6], marker="s", c="C1")
#
#     img.set_label(0,0,0,"$g_1$")
#     img.set_label(0,0,1,"$\hat g_1$")
#     img.set_label(0,1,0,"$g_2$")
#     img.set_label(0,1,1,"$\hat g_2$")
#     img.axs[0][0].legend(fontsize=img.legend_size)
#     img.axs[0][1].legend(fontsize=img.legend_size)
#     img.save_img("Galaxies_check.png")
#     img.close_img()
#
#     # sigma ratio
#     img = Image_Plot()
#     img.subplots(1,1)
#     img.axs[0][0].plot(range(result_col), result[1]/result[5], c="C2",label="$\sigma$ ratio of $g_1$")
#     img.axs[0][0].plot(range(result_col), result[3]/result[7], c="C1",label="$\sigma$ ratio of $g_2$")
#     img.set_label(0,0,0,"$\sigma_{ori} / \sigma_{new}$ ratio")
#     img.axs[0][0].legend(fontsize=img.legend_size)
#     img.save_img("Galaxies_check_sigma.png")
#     img.close_img()