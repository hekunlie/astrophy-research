import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
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
    # print(num, n1, n2)
    return numpy.sum(xi[:len(xi) - ig_num]) * 0.5, xi

def set_bin(data, bin_num, bound_scale):
    temp_data = numpy.sort(numpy.abs(data))
    bin_size = len(temp_data)/bin_num*2
    bins = numpy.array([temp_data[int(i*bin_size)] for i in range(1, int(bin_num / 2))])
    bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
    bound = numpy.max(numpy.abs(data)) * bound_scale
    bins = numpy.append(-bound, numpy.append(bins, bound))
    return bins

fmt='%2.f%%'
fig_x = 9
fig_y = fig_x*4/6
fonts = 20
xy_lb_size = 22
xy_tick_size = xy_lb_size - 5
legend_size =xy_tick_size
axis_linewidth = 2
plt_line_width = 2
cap_size = 5
tick_len = 6
figs = (fig_x*4, fig_y*2)
fig = plt.figure(figsize=figs)
axs = []
for i in range(8):
    ax = fig.add_subplot(241 + i)
    axs.append(ax)

sources = [int(argv[2]), int(argv[3])]
bin_num = 12

total_path = "/mnt/ddnfs/data_users/hkli/selection_bias_64/"

shear_cata = total_path + "parameters/shear.npz"
shear = numpy.load(shear_cata)


locs = ["upper left", "lower left"]
for m in range(2):
    source = sources[m]
    g1 = shear["arr_0"][source]
    g2 = shear["arr_1"][source]

    # point source
    fourier_path = total_path + "result/data/data_%d.hdf5"%source
    f_h5 = h5py.File(fourier_path, "r")
    es_data = f_h5["/data"].value
    f_h5.close()
    if argv[1] == "e1":
        idx = 0
    else:
        idx = 1
    e1 = es_data[:,idx]
    # ch = numpy.random.choice(range(len(e1)), 10000, replace=False)

    # galsim
    # para_path = total_path + "parameters/para_%d.hdf5"%source
    # para_h5 = h5py.File(para_path, "r")
    # mag_t = para_h5["/mag"].value
    # e1 = para_h5["/e1"].value
    # # e1 = para_h5["/%s"%argv[2]].value
    # para_h5.close()

    bins = tool_box.set_bin(e1, bin_num, 1.02)

    fourier_path = total_path + "result/data/data_1.5sig/data_%d.hdf5"%source
    f_h5 = h5py.File(fourier_path, "r")
    f_data = f_h5["/data"].value
    f_h5.close()

    sex_path = total_path + "result/data/sex2_1.5/sex_%d.npz"%source
    s_data = numpy.load(sex_path)["arr_0"]

    detect_f = f_data[:, -1] > -1
    detect_s = s_data[:, 0] > 0

    # P(k=0)
    pk0 = f_data[:, 4]/64/60
    print(pk0.max(), pk0.min())

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

    cut_data = [pk0, snr, snr_auto, -mag]

    labels = ["P$_{k0}$", "SNR$_S$", "SNR$_A$", "MAG_AUTO", "MAG_TRUE_f", "MAG_TRUE_s"]
    cuts_num = 10

    plt_ys = [0, 0]
    plt_xs = [0, 0]

    for i in range(4):
        ax_tag = m*4+i
        xi_st, xi_t = G_bin(e1, bins)
        axs[ax_tag].errorbar(range(1, len(xi_t)+1),xi_t, marker="s", label="Total")
        if i == 0 or i == 4:
            idx = detect_f
            # ax.hist(e1[detect_f], bins, alpha=0.8, label="Fourier Quad")
            xi_s, xi = G_bin(e1[detect_f], bins)
            axs[ax_tag].errorbar(range(1, len(xi)+1), xi, marker="s", label="Detected")
        else:
            idx = detect_s
            # ax.hist(e1[detect_s], bins, alpha=0.8, label="SExtractor")
            xi_s, xi = G_bin(e1[detect_s], bins)
            axs[ax_tag].errorbar(range(1,len(xi)+1), xi,  marker="s",label="Detected")
        cut_data_sort = numpy.sort(cut_data[i][idx])
        cut_data_step = int(len(cut_data_sort) / cuts_num)
        cutoff_scalar = [cut_data_sort[i * cut_data_step] for i in range(cuts_num)]
        for j in range(2):
            # continue
            cut_s = cutoff_scalar[3+j*4]
            if j == 0:
                if i == 3:
                    lb = labels[i] + "$\leq$ %.2f (70%% remain)" %(-cut_s)
                else:
                    lb = labels[i] + "$\geq$ %.2f (70%% remain)"%cut_s
            else:
                if i == 3:
                    lb = labels[i] + "$\leq$ %.2f (30%% remain)"%(-cut_s)
                else:
                    lb = labels[i] + "$\geq$ %.2f (30%% remain)"%cut_s

            idx_scale = cut_data[i] >= cut_s
            xi_s, xi = G_bin(e1[idx&idx_scale], bins)
            axs[ax_tag].errorbar(range(1,len(xi)+1), xi,  marker="s", label=lb)
        ys = axs[ax_tag].set_ylim()
        xs = axs[ax_tag].set_xlim()
        if ys[1] > plt_ys[1]:
            plt_ys[1] = ys[1]
        if ys[0] < plt_ys[0]:
            plt_ys[0] = ys[0]
        if xs[1] > plt_xs[1]:
            plt_xs[1] = xs[1]
        if xs[0] < plt_xs[0]:
            plt_xs[0] = xs[0]
        axs[ax_tag].legend(loc=locs[m], fontsize=legend_size)
    # ax.set_ylim(ys[0], ys[1])
    # print(plt_xs, plt_ys)

    for i in range(4):
        ax_tag = m * 4 + i
        axs[ax_tag].plot([plt_xs[0],plt_xs[1]], [0,0], linestyle="--", c="grey")
        axs[ax_tag].set_ylim(plt_ys[0]-3, plt_ys[1]+2)
        if i > 0:
            axs[ax_tag].set_yticklabels([])
        else:
            if "e1" == argv[1]:
                e_label = "$\chi^2$ of e$_1$ (g1=%.2f)"%g1
            else:
                e_label = "$\chi^2$ of e$_2$ (g2=%.2f)"%g2
            axs[ax_tag].set_ylabel(e_label, fontsize=xy_lb_size)
        if m == 1:
            axs[ax_tag].set_xlabel("Bin label", fontsize=xy_lb_size)
        else:
            axs[ax_tag].set_xticklabels([])
        axs[ax_tag].tick_params(direction='in', labelsize=xy_tick_size, top=True, right=True,pad=5)
        for axis in ["bottom", "left", "top", "right"]:
            # the line width of the frame
            axs[ax_tag].spines[axis].set_linewidth(axis_linewidth)
        axs[ax_tag].xaxis.set_tick_params(which="both",direction="in",length=tick_len, width=axis_linewidth)
        axs[ax_tag].yaxis.set_tick_params(which="major",direction="in",length=tick_len, width=axis_linewidth)
        axs[ax_tag].yaxis.set_tick_params(which="minor",direction="in",length=tick_len-2, width=axis_linewidth)
    print(g1, g2)
plt.subplots_adjust(wspace=0, hspace=0)
plt.savefig("/home/hkli/work/selection_bias/sym_mc_plot/pics/%s_%s_%s.pdf"%(argv[1], argv[2], argv[3]),bbox_inches='tight')

