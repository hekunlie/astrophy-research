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


cmd = "e_1"
# sources = [2, 8]
sources = [2, 5]
bin_num = 12
cuts_num = 10

# total_path = "/mnt/ddnfs/data_users/hkli/selection_bias_64_1/"
# total_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/galsim_dimmer/"
total_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/pts_dimmer/"

pic_nm = "chisq_pts.pdf"
pic_nm_png = "chisq_pts.png"

sex_filter_name = "sex2_1.5"

shear_cata = total_path + "parameters/shear.npz"
shear = numpy.load(shear_cata)

fmt = '%2.f%%'

img = Image_Plot(fig_x=5, fig_y=4, xpad=0.45,ypad=0.22)
img.subplots(2, 2)
img.axis_type(0, "major", tick_len=8, tick_width=2)
img.axis_type(1, "major", tick_len=8, tick_width=2)

locs = ["upper left", "upper left","upper left","upper left"]
legend_loc = [(0.02, 0.95), (0.02, 0.95), (0.02, 0.95), (0.02, 0.66)]
legend_loc_share_ax = [(0.02, 0.62), (0.02, 0.62), (0.02, 0.62), (0.02, 0.35)]
ylims = [(-9, 85),(-24, 285),(-9, 85),(-850, 85)]
for m in range(1):
    source = sources[m]
    g1 = shear["arr_0"][source]
    g2 = shear["arr_1"][source]

    # point source
    fourier_path = total_path + "result/data/data_%d.hdf5" % source
    f_h5 = h5py.File(fourier_path, "r")
    es_data = f_h5["/data"].value
    f_h5.close()
    if cmd == "e_1":
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

    bins = tool_box.set_bin(e1, bin_num, 1.01)
    num_ini = numpy.histogram(e1, bins)[0]

    cut_data = []

    selection_labels = ["P$_{k0}$","SNR$_S$","MAG_AUTO","Resolution"]
    data_paths = [total_path + "result/data/%s/flux2_ex1_%d.hdf5" % (sex_filter_name, source),
                  total_path + "result/data/%s/snr_sex_%d.hdf5" % (sex_filter_name, source),
                  total_path + "result/data/%s/mag_auto_%d.hdf5" % (sex_filter_name, source),
                  total_path + "result/data/%s/rfactor_%d.hdf5" % (sex_filter_name, source)]

    for data_path in data_paths:
        h5f = h5py.File(data_path, "r")
        data = h5f["/data"].value
        h5f.close()
        if len(data.shape) == 2:
            data.shape = (data.shape[0],)
        if "flux2_" in data_path:
            data = data/64/60
        cut_data.append(data)

    mask_path = total_path + "result/data/%s/mask_%d.hdf5" % (sex_filter_name, source)
    h5f = h5py.File(mask_path, "r")
    mask = h5f["/data"].value
    h5f.close()
    # detected
    idx_detect = mask > 0
    print("Detecton: ",idx_detect.sum())
    plt_ys = [0, 0]
    plt_xs = [0, 0]

    # the chi squared of total sample
    chi_total_sum, chi_total, num_total = G_bin(e1, bins)
    # the chi squared of detected sample
    chi_detect_sum, chi_detect, num_detect = G_bin(e1[idx_detect], bins)

    ellip_scale = (bins[1:] + bins[:bin_num])[int(bin_num/2):]/2
    dscale = (ellip_scale[-1] - ellip_scale[0])*0.05
    xscale_lim = (ellip_scale[0] - ellip_scale[0]*0.2, ellip_scale[-1] + ellip_scale[-1]*0.1)
    # print(ellip_scale)
    # print(xscale_lim)

    for i in range(len(cut_data)):
        ny,nx = divmod(i,2)
        # the chi squared of total sample
        img.axs[ny][nx].scatter(range(1, len(chi_total) + 1), chi_total, marker="s", label="Total", facecolor="none",
                              linewidths=img.plt_line_width, color="C0", s=80)
        img.axs[ny][nx].plot(range(1, len(chi_total) + 1), chi_total, color="C0")

        # the chi squared of detected sample
        img.axs[ny][nx].scatter(range(1, len(chi_detect) + 1), chi_detect, marker="s", label="Detected", facecolor="none",
                              linewidths=img.plt_line_width,   color="C1", s=80)
        img.axs[ny][nx].plot(range(1, len(chi_detect) + 1), chi_detect, color="C1")

        # plot the mean ellipticity of each bin
        share_ax = img.share_axis(ny, nx, 1)
        ylabel = "|$%s$|" % cmd
        share_ax.scatter(range(1, len(chi_detect) + 1), ellip_scale, color="C4",marker="o",facecolor="none",
                         linewidth=img.plt_line_width,label=ylabel, s=80)
        share_ax.plot(range(1, len(chi_detect) + 1), ellip_scale, color="C4",linewidth=img.plt_line_width)
        share_ax.set_ylabel(ylabel, fontsize=img.xy_lb_size)
        share_ax.legend(loc=locs[i], fontsize=img.legend_size, frameon=False,bbox_to_anchor=legend_loc_share_ax[i])


        # img.axs[m][i].scatter(ellip_scale, xi, marker="s", label="Detected", facecolor="none",
        #                       linewidths=img.plt_line_width,   color="C1", s=80)
        # img.axs[m][i].plot(ellip_scale, xi, color="C1")
        #
        # img.axs[m][i].scatter(range(1, len(num_s) + 1), num_s/num_ini, marker="s", label="Detected", facecolor="none",
        #                       linewidths=img.plt_line_width,   color="C1", s=80)
        # img.axs[m][i].plot(range(1, len(num_s) + 1), num_s/num_ini, color="C1")

        # plot the chi squared of sample with 30% and 70% source cutoff
        cut_data_sort = numpy.sort(cut_data[i][idx_detect])
        cut_data_step = int(len(cut_data_sort) / cuts_num)
        cutoff_scalar = [cut_data_sort[i * cut_data_step] for i in range(cuts_num)]
        print(selection_labels[i])
        print(cutoff_scalar)
        for j in range(2):
            # continue
            cut_s = cutoff_scalar[3 + j * 4]

            if "MAG" in selection_labels[i]:
                lb = selection_labels[i] + "$\leq %.2f$" % (-cut_s)
            else:
                lb = selection_labels[i] + "$\geq %.2f$" % cut_s

            idx_scale = cut_data[i] >= cut_s

            chi_cut_sum, chi_cut, num_cut = G_bin(e1[idx_detect & idx_scale], bins)

            img.axs[ny][nx].scatter(range(1, len(chi_cut) + 1), chi_cut, marker="s", label=lb, facecolor="none",
                                  linewidths=img.plt_line_width, color="C%d" % (j + 2), s=80)
            img.axs[ny][nx].plot(range(1, len(chi_cut) + 1), chi_cut, color="C%d" % (j + 2))

            # img.axs[m][i].scatter(ellip_scale, xi, marker="s", label=lb, facecolor="none",
            #                       linewidths=img.plt_line_width, color="C%d" % (j + 2), s=80)
            # img.axs[m][i].plot(ellip_scale, xi, color="C%d" % (j + 2))

            # img.axs[m][i].scatter(range(1, len(num_s) + 1), num_s/num_ini, marker="s", label=lb, facecolor="none",
            #                       linewidths=img.plt_line_width, color="C%d" % (j + 2), s=80)
            # img.axs[m][i].plot(range(1, len(num_s) + 1), num_s/num_ini, color="C%d" % (j + 2))

        img.axs[ny][nx].set_ylim(ylims[i])
        img.axs[ny][nx].set_xlim()

        img.axs[ny][nx].legend(loc=locs[i], fontsize=img.legend_size, frameon=False, bbox_to_anchor=legend_loc[i])


        e_label = "S$_{%s}$, $g_1$=%.2f"%(cmd, g1)
        img.axs[ny][nx].set_ylabel(e_label, fontsize=img.xy_lb_size)

        dy = (plt_ys[1] - plt_ys[0])*0.05
        # img.axs[m][i].set_ylim(plt_ys[0] - dy, plt_ys[1] + dy)
        img.axs[ny][nx].set_xlabel("Bin label", fontsize=img.xy_lb_size)

    print(g1, g2)
# img.figure.subplots_adjust(wspace=0, hspace=0)
img.save_img(pic_nm)
img.save_img(pic_nm_png)




