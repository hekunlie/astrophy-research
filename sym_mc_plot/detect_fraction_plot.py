import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import ScalarFormatter
import tool_box
import h5py
from Fourier_Quad import Fourier_Quad
from sys import argv

def G_bin(G_h, bins, ig_num=0):  # checked 2017-7-9!!!
    r"""
    to calculate the symmetry the shear estimators
    :param g: estimators from Fourier quad, 1-D numpy array
    :param nu: N + U for g1, N - U for g2, 1-D numpy array
    :param g_h: pseudo shear (guess)
    :param bins: bin of g for calculation of the symmetry, 1-D numpy array
    :param ig_num: the number of inner grid of bin to be neglected
    :return: chi square
    """
    bin_num = len(bins) - 1
    inverse = range(int(bin_num / 2 - 1), -1, -1)
    num = numpy.histogram(G_h, bins)[0]
    n1 = num[0:int(bin_num / 2)][inverse]
    n2 = num[int(bin_num / 2):]
    xi = (n1 - n2) ** 2 / (n1 + n2)
    delta = n1 - n2
    idx = delta < 0
    xi[idx] = -xi[idx]
    print(num, n1, n2)
    return numpy.sum(xi[:len(xi) - ig_num]) * 0.5, xi

def set_bin(data, bin_num, bound_scale):
    temp_data = numpy.sort(numpy.sort(data))
    bin_size = len(temp_data)/bin_num*2
    bins = numpy.array([temp_data[int(i*bin_size)] for i in range(1, int(bin_num / 2))])
    bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
    bound = numpy.max(numpy.abs(data)) * bound_scale
    bins = numpy.append(-bound, numpy.append(bins, bound))
    return bins

fmt='%2.f%%'
fig_x = 6
fig_y = fig_x*4/6
figs = (fig_x*4, fig_y)
fonts = 20
xy_lb_size = 18
legend_size = fonts - 4
axis_linewidth = 1.2
plt_line_width = 2
cap_size = 5
# xticks = mtick.FormatStrFormatter(fmt)

fig = plt.figure(figsize=figs)

source = int(argv[1])

total_path = "/mnt/ddnfs/data_users/hkli/simu_test/"

shear_cata = total_path + "parameters/shear.npz"
shear = numpy.load(shear_cata)
g1 = shear["arr_0"][source]
g2 = shear["arr_1"][source]


para_path = total_path + "parameters/para_%d.hdf5"%source
para_h5 = h5py.File(para_path, "r")
mag_t = para_h5["/mag"].value
# e1 = para_h5["/e1"].value
e1 = para_h5["/%s"%argv[2]].value
para_h5.close()
bins = tool_box.set_bin(e1, int(argv[3]), 1.02)
# bins = numpy.linspace(-0.81,0.81, int(argv[1]))
print(bins)
fourier_path = total_path + "result/data/data_1.5sig/data_%d.hdf5"%source
f_h5 = h5py.File(fourier_path, "r")
f_data = f_h5["/data"].value
f_h5.close()

sex_path = total_path + "result/data/sex2_1.5/sex_%d.npz"%source
s_data = numpy.load(sex_path)["arr_0"]

detect_f = f_data[:, -1] > -1
detect_s = s_data[:, 0] > 0

# P(k=0)
pk0 = f_data[:, 4]

# snr_auto
flux_auto = s_data[:, 1]
flux_err = s_data[:, 2]
idx_s = flux_err == 0
flux_err[idx_s] = 1
snr_auto = flux_auto / flux_err

# snr
snr = s_data[:, 0]

# mag_auto
mag = s_data[:, 3]


cut_data = [pk0, snr_auto, snr, -mag, -mag_t, -mag_t]

labels = ["P(k=0)", "SNR_AUTO", "SNR", "MAG_AUTO", "MAG_TRUE_f", "MAG_TRUE_s"]
cuts_num = 10
# rng = numpy.random.RandomState(123)
#
# ch = rng.choice(range(4000000), 30000, replace=False)
axs = []
plt_ys = [0,0]
plt_xs = [0,0]
for i in range(4):
    ax = fig.add_subplot(141 + i)
    ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
    for axis in ["bottom", "left", "top", "right"]:
        # the line width of the frame
        ax.spines[axis].set_linewidth(axis_linewidth)
    ax.xaxis.set_tick_params(which="both", direction="in", length=6, width=axis_linewidth)
    axs.append(ax)

for i in range(4):

    xi_st, xi_t = G_bin(e1, bins)
    axs[i].errorbar(range(len(xi_t)),xi_t, marker="s", label="TOTAL SAMPLE")
    if i == 0 or i == 4:
        idx = detect_f
        # ax.hist(e1[detect_f], bins, alpha=0.8, label="Fourier Quad")
        xi_s, xi = G_bin(e1[detect_f], bins)
        axs[i].errorbar(range(len(xi)), xi, marker="s", label="Fourier Quad")
    else:
        idx = detect_s
        # ax.hist(e1[detect_s], bins, alpha=0.8, label="SExtractor")
        xi_s, xi = G_bin(e1[detect_s], bins)
        axs[i].errorbar(range(len(xi)), xi,  marker="s",label="SExtractor")


    # ax.scatter(e1[idx][ch], cut_data[i][idx][ch], s=0.05)
    # ax.set_xlabel("e2")
    # xs = ax.set_xlim()
    # ys = ax.set_ylim()
    # if i < 4:
    #     y_span = ys[0] + (cutoff_scalar[6] - ys[0])*1.3
    #     for j in range(1,4):
    #         cut_line = cutoff_scalar[j*2]
    #         ax.plot([xs[0], xs[1]],[cut_line, cut_line],linestyle="--",label="%d0%%"%(j*2))
    #     ax.legend(ncol=4)
        # ax.set_ylim(ys[0],y_span)
    # print(labels[i], cutoff_scalar)

    cut_data_sort = numpy.sort(cut_data[i][idx])
    cut_data_step = int(len(cut_data_sort) / cuts_num)
    cutoff_scalar = [cut_data_sort[i * cut_data_step] for i in range(cuts_num)]
    for j in range(2):
        # continue
        if j == 0:
            lines = "-"
            lb = labels[i] + " 30%"
        else:
            lines = "--"
            lb = labels[i] + " 60%"
        cut_s = cutoff_scalar[3+j*3]
        idx_scale = cut_data[i] >= cut_s
        xi_s, xi = G_bin(e1[idx&idx_scale], bins)
        print(lb, xi_s)
        axs[i].errorbar(range(len(xi)), xi,  marker="s", label=lb)
        # tag = numpy.where(cut_s == cut_data[i][idx])
        # mag_t_ch = mag_t[idx][:,0][tag[0].tolist()]
        # mag_t_max, mag_t_min = mag_t_ch.max(), mag_t_ch.min()
        # mag_refer = mag_t[idx][tag]

        # print(labels[i], cut_s, tag)
        # print(labels[i], mag_t_max, mag_t_min)

        # if mag_t_max == mag_t_min:
        #     mag_refer = mag_t_max
        #     ax.plot([mag_refer, mag_refer], [ys[0], ys[1]],  linestyle=lines, label=lb)
        # else:
        #     ax.plot([mag_t_min, mag_t_min], [ys[0], ys[1]],  linestyle=lines, label=lb+"_min")
        #     ax.plot([mag_t_max, mag_t_max], [ys[0], ys[1]],  linestyle=lines, label=lb+"_max")
    ys = axs[i].set_ylim()
    xs = axs[i].set_xlim()
    if ys[1] > plt_ys[1]:
        plt_ys[1] = ys[1]
    if ys[0] < plt_ys[0]:
        plt_ys[0] = ys[0]
    if xs[1] > plt_xs[1]:
        plt_xs[1] = xs[1]
    if xs[0] < plt_xs[0]:
        plt_xs[0] = xs[0]
    axs[i].legend()
# ax.set_ylim(ys[0], ys[1])
for i in range(4):
    axs[i].plot([plt_xs[0],plt_xs[1]], [0,0], linestyle="--", c="grey")
    axs[i].set_ylim(plt_ys[0], plt_ys[1])
    # if i > 0:
    #     axs[i].set_yticks([])
plt.suptitle("g1:%1.4f, g2:%1.4f"%(g1, g2))
plt.subplots_adjust(wspace=0)
plt.savefig("/home/hkli/work/selection_bias/sym_mc_plot/pics/%s_%s_%s_chisq.png"%(argv[1],argv[2], argv[3]),bbox_inches='tight')
