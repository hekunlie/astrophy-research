import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
path.append('D:/GitHub/astrophy-research/mylib')
import numpy
import matplotlib.ticker as mtick
from plot_tool import Image_Plot
import h5py

#  the m & c of different filters from point source


matplotlib.rcParams["font.family"] = "serif"

source_f = "pts_bright"
# final pic name
pic_nm = "pts_filters_f_m1.pdf"

# the faint source
data_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/cuts/sym/" % source_f

sex_filters = ["sex2_1.5", "sex4_1.5", "sex2_4", "sex4_4"]
names = ["P$_{k0}$", "MAG_AUTO", "MAG$_{true}$", "SNR$_S$", "SNR$_A$"]
files = "mag_auto"
ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]
ylabels = ["$m_1 \\times 10^2$", "$m_2 \\times 10^2$", "$m_1 \\times 10^2$", "$m_2 \\times 10^2$"]


x_tick = [i * cuts_num for i in range(0,ch_num,2)]
fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

img = Image_Plot(fig_x=9,fig_y=6)
img.subplots(1,1)

dys = (-0.1, 0.46)
legend_pos = (0.02, 0.96)

col = 0

gauss_size = [2, 4, 2, 4]
sig_scale = [1.5, 1.5, 4, 4]
color = ["C0", "C4", "C1", "C3"]
for tag, filter_name in enumerate(sex_filters):
    print(filter_name, files)
    if tag < 2:
        data = numpy.load(data_path + "%s/%s/total.npz" % (filter_name, files))
        lb = "$MAG\_AUTO, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = data['arr_%d' % col][:, ch]
        img.axs[0][0].errorbar(x_coord, 100 * (mc[0] - 1), 100 * mc[1], c=color[tag], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=lb, marker="s", fillstyle="none")

        data = numpy.load(data_path + "%s/flux2_ex1/total.npz" % filter_name)
        lb = "$P_{k0}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = data['arr_%d' % col][:, ch]
        img.axs[0][0].errorbar(x_coord, 100 * (mc[0] - 1), 100 * mc[1], c=color[tag + 2], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=lb, marker="p")
    else:
        data = numpy.load(data_path + "%s/%s/total.npz" % (filter_name, files))
        lb = "$MAG\_AUTO, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = data['arr_%d' % col][:, ch]
        img.axs[0][0].errorbar(x_coord, 100 * (mc[0] - 1), 100 * mc[1], c=color[tag - 2], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=lb, marker="s", linestyle="--", alpha=0.65, fillstyle="none")

        data = numpy.load(data_path + "%s/flux2_ex1/total.npz" % filter_name)
        lb = "$P_{k0}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = data['arr_%d' % col][:, ch]
        img.axs[0][0].errorbar(x_coord, 100 * (mc[0] - 1), 100 * mc[1], c=color[tag], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=lb, marker="p", linestyle="--", alpha=0.65)

xs = img.axs[0][0].set_xlim()
ys = img.axs[0][0].set_ylim()
img.axs[0][0].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width, c="grey", linestyle="--")
img.axs[0][0].set_xlim(xs[0], xs[1])
img.axs[0][0].set_ylim(ys[0] + dys[0], ys[1] + dys[1])
img.axs[0][0].xaxis.set_major_formatter(xticks)
img.axs[0][0].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
img.axs[0][0].set_ylabel(ylabels[col], fontsize=img.xy_lb_size)
img.axs[0][0].set_xticks(x_tick)
img.axs[0][0].legend(ncol=2, loc="upper left", fontsize=img.legend_size-2, frameon=False)
img.save_img(pic_nm)
