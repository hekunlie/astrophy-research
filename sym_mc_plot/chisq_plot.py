import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/' % my_home)
import numpy
import matplotlib.pyplot as plt
import tool_box
import h5py
import matplotlib
from plot_tool import Image_Plot



matplotlib.style.use('default')
plt.rcParams['font.family'] = 'serif'


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
    # xi = (n1 - n2) / (n1 + n2)#*numpy.abs(n1-n2)
    xi = (n1 - n2) ** 2 / (n1 + n2)
    delta = n1 - n2
    idx = delta < 0
    xi[idx] = -xi[idx]
    # print(num, n1, n2)
    return numpy.sum(xi[:len(xi) - ig_num]) * 0.5, xi, num


def set_bin(data, bin_num, bound_scale):
    temp_data = numpy.sort(numpy.abs(data))
    bin_size = len(temp_data) / bin_num * 2
    bins = numpy.array([temp_data[int(i * bin_size)] for i in range(1, int(bin_num / 2))])
    bins = numpy.sort(numpy.append(numpy.append(-bins, [0.]), bins))
    bound = numpy.max(numpy.abs(data)) * bound_scale
    bins = numpy.append(-bound, numpy.append(bins, bound))
    return bins


cmd = "e1"
# sources = [2, 8]
sources = [2, 5]
bin_num = 12

# total_path = "/mnt/ddnfs/data_users/hkli/selection_bias_64_1/"
# total_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/galsim_dimmer/"
total_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/pts_dimmer/"

pic_nm = "chisq_pts.pdf"

shear_cata = total_path + "parameters/shear.npz"
shear = numpy.load(shear_cata)

fmt = '%2.f%%'
fig_x = 6
fig_y = fig_x * 4 / 6
img = Image_Plot(fig_x=fig_x, fig_y=fig_y)
img.subplots(2,4)


locs = ["upper left", "lower left"]
legend_loc = [(0.02, 0.97), (0.02, 0.08)]
for m in range(2):
    source = sources[m]
    g1 = shear["arr_0"][source]
    g2 = shear["arr_1"][source]

    # point source
    fourier_path = total_path + "result/data/data_%d.hdf5" % source
    f_h5 = h5py.File(fourier_path, "r")
    es_data = f_h5["/data"].value
    f_h5.close()
    if cmd == "e1":
        idx = 0
    else:
        idx = 1
    e1 = es_data[:, idx]

    # # galsim
    # para_path = total_path + "parameters/para_%d.hdf5"%source
    # para_h5 = h5py.File(para_path, "r")
    # mag_t = para_h5["/mag"].value
    # e1 = para_h5["/e1"].value
    # para_h5.close()

    bins = tool_box.set_bin(e1, bin_num, 1.02)
    num_ini = numpy.histogram(e1, bins)[0]

    fourier_path = total_path + "result/data/data_1.5sig/data_%d.hdf5" % source
    f_h5 = h5py.File(fourier_path, "r")
    f_data = f_h5["/data"].value
    f_h5.close()

    sex_path = total_path + "result/data/sex2_1.5/sex_%d.npz" % source
    s_data = numpy.load(sex_path)["arr_0"]

    detect_f = f_data[:, -1] > -1
    detect_s = s_data[:, 0] > 0

    # P(k=0)
    pk0 = f_data[:, 4]/64/60

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

    labels = ["P$_{k0}$", "SNR$_S$", "SNR$_A$", "MAG_AUTO"]
    cuts_num = 10

    plt_ys = [0, 0]
    plt_xs = [0, 0]

    xi_st, xi_t, num_t = G_bin(e1, bins)
    ellip_scale = (bins[1:] + bins[:bin_num])[int(bin_num/2):]/2

    dscale = (ellip_scale[-1] - ellip_scale[0])*0.05
    xscale_lim = (ellip_scale[0] - ellip_scale[0]*0.2, ellip_scale[-1] + ellip_scale[-1]*0.1)
    print(ellip_scale)
    print(xscale_lim)
    for i in range(4):
        # axs[ax_tag].errorbar(range(1, len(xi_t)+1),xi_t, marker="s", label="Total")
        img.axs[m][i].scatter(range(1, len(xi_t) + 1), xi_t, marker="s", label="Total", facecolor="none",
                              linewidths=img.plt_line_width, color="C0", s=80)
        img.axs[m][i].plot(range(1, len(xi_t) + 1), xi_t, color="C0")

        # img.axs[m][i].scatter(ellip_scale, xi_t, marker="s", label="Total", facecolor="none",
        #                       linewidths=img.plt_line_width, color="C0", s=80)
        # img.axs[m][i].plot(ellip_scale, xi_t, color="C0")

        # img.axs[m][i].scatter(range(1, len(num_t) + 1), num_t/num_ini, marker="s", label="Total", facecolor="none",
        #                       linewidths=img.plt_line_width, color="C0", s=80)
        # img.axs[m][i].plot(range(1, len(num_t) + 1), num_t/num_ini, color="C0")

        if i == 0:
            idx = detect_s
            # ax.hist(e1[detect_f], bins, alpha=0.8, label="Fourier Quad")

            xi_s, xi, num_s = G_bin(e1[detect_s], bins)
            #             axs[ax_tag].errorbar(range(1, len(xi)+1), xi, marker="s", label="Detected")
        else:
            idx = detect_s
            # ax.hist(e1[detect_s], bins, alpha=0.8, label="SExtractor")

            xi_s, xi, num_s = G_bin(e1[detect_s], bins)
            #             axs[ax_tag].errorbar(range(1,len(xi)+1), xi,  marker="s",label="Detected")

        img.axs[m][i].scatter(range(1, len(xi) + 1), xi, marker="s", label="Detected", facecolor="none",
                              linewidths=img.plt_line_width,   color="C1", s=80)
        img.axs[m][i].plot(range(1, len(xi) + 1), xi, color="C1")

        # img.axs[m][i].scatter(ellip_scale, xi, marker="s", label="Detected", facecolor="none",
        #                       linewidths=img.plt_line_width,   color="C1", s=80)
        # img.axs[m][i].plot(ellip_scale, xi, color="C1")
        #
        # img.axs[m][i].scatter(range(1, len(num_s) + 1), num_s/num_ini, marker="s", label="Detected", facecolor="none",
        #                       linewidths=img.plt_line_width,   color="C1", s=80)
        # img.axs[m][i].plot(range(1, len(num_s) + 1), num_s/num_ini, color="C1")

        # plot the chi squared of cutoff
        cut_data_sort = numpy.sort(cut_data[i][idx])
        cut_data_step = int(len(cut_data_sort) / cuts_num)
        cutoff_scalar = [cut_data_sort[i * cut_data_step] for i in range(cuts_num)]
        for j in range(2):
            # continue
            cut_s = cutoff_scalar[3 + j * 4]
            if j == 0:
                if i == 3:
                    lb = labels[i] + "$\leq %.2f$" % (-cut_s)
                else:
                    lb = labels[i] + "$\geq %.2f$" % cut_s
            else:
                if i == 3:
                    lb = labels[i] + "$\leq %.2f$" % (-cut_s)
                else:
                    lb = labels[i] + "$\geq %.2f$" % cut_s

            idx_scale = cut_data[i] >= cut_s

            xi_s, xi, num_s = G_bin(e1[idx & idx_scale], bins)
            # axs[ax_tag].errorbar(range(1,len(xi)+1), xi,  marker="s", label=lb,facecolor="none",linewidths=img.plt_line_width)
            img.axs[m][i].scatter(range(1, len(xi) + 1), xi, marker="s", label=lb, facecolor="none",
                                  linewidths=img.plt_line_width, color="C%d" % (j + 2), s=80)
            img.axs[m][i].plot(range(1, len(xi) + 1), xi, color="C%d" % (j + 2))

            # img.axs[m][i].scatter(ellip_scale, xi, marker="s", label=lb, facecolor="none",
            #                       linewidths=img.plt_line_width, color="C%d" % (j + 2), s=80)
            # img.axs[m][i].plot(ellip_scale, xi, color="C%d" % (j + 2))

            # img.axs[m][i].scatter(range(1, len(num_s) + 1), num_s/num_ini, marker="s", label=lb, facecolor="none",
            #                       linewidths=img.plt_line_width, color="C%d" % (j + 2), s=80)
            # img.axs[m][i].plot(range(1, len(num_s) + 1), num_s/num_ini, color="C%d" % (j + 2))

        ys = img.axs[m][i].set_ylim()
        xs = img.axs[m][i].set_xlim()
        if ys[1] > plt_ys[1]:
            plt_ys[1] = ys[1]
        if ys[0] < plt_ys[0]:
            plt_ys[0] = ys[0]
        if xs[1] > plt_xs[1]:
            plt_xs[1] = xs[1]
        if xs[0] < plt_xs[0]:
            plt_xs[0] = xs[0]
        img.axs[m][i].legend(loc=locs[m], fontsize=img.legend_size, frameon=False, bbox_to_anchor=legend_loc[m])
    # ax.set_ylim(ys[0], ys[1])
    # print(plt_xs, plt_ys)

    for i in range(4):
        # img.axs[m][i].plot([plt_xs[0], plt_xs[1]], [0, 0], linestyle="--", c="grey")


        if i > 0:
            img.axs[m][i].set_yticklabels([])
        else:
            # if m == 0:
            #     axs[ax_tag].set_yticks(numpy.arange(0, 300, 50))
            # else:
            #     axs[ax_tag].set_yticks(numpy.arange(-250, 50, 50))
            if "e1" == cmd:
                e_label = "S$_{e_1}$, $g_1$=%.2f"%g1
            else:
                e_label = "S($e_2$)"
            img.axs[m][i].set_ylabel(e_label, fontsize=img.xy_lb_size)
        dy = (plt_ys[1] - plt_ys[0])*0.05
        img.axs[m][i].set_ylim(plt_ys[0] - dy, plt_ys[1] + dy)
        if m == 1:
            text_y = -100
            img.axs[m][i].set_xlabel("Bin label", fontsize=img.xy_lb_size)
        else:
            # img.axs[m][i].set_ylim(plt_ys[0], plt_ys[1])
            text_y = 100
            img.axs[m][i].set_xticklabels([])

        # img.axs[m][i].text(1.2, text_y, '$g_1$=%.2f'%g1, color='black', ha='left', va='center',  fontsize=img.legend_size)

        # img.axs[m][i].set_xscale("log")
        # img.axs[m][i].set_xlim(xscale_lim[0], xscale_lim[1])
    print(g1, g2)
img.figure.subplots_adjust(wspace=0, hspace=0)
img.save_img(pic_nm)




