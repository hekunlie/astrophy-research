import matplotlib
# matplotlib.use("Agg")
import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
# path.append('%s/work/mylib/'%my_home)
path.append('D:/GitHub/astrophy-research/mylib')
import numpy
import matplotlib.ticker as mtick
from plot_tool import Image_Plot
import h5py

#  the m & c of different filters from point source


matplotlib.rcParams["font.family"] = "serif"


# final pic name
pic_nm = "E:/comparison.pdf"

# the faint source
# data_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/cuts/sym/" % source_f
data_path = ["D:/cut/large_g/galsim_dimmer/result/cuts_pi_all_sample_w_tflux_sq/sym/",
             "D:/cut/large_g/pts_dimmer/result/cuts_pi_all_sample_w_tflux_sq/sym/"]


sex_filters = ["sex2_1.5", "sex4_1.5", "sex2_4", "sex4_4"]

files = "mag_auto"
ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]
ylabels = ["$m_1 \\times 10^2$", "$m_2 \\times 10^2$", "$m_1 \\times 10^2$", "$m_2 \\times 10^2$"]


x_tick = [i * cuts_num for i in range(0,ch_num,2)]
fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

img = Image_Plot(fig_x=6,fig_y=4,legend_size=12, plt_line_width=2,xpad=0.20)
img.subplots(1,2)
img.axis_type(0,"major",tick_len=6, tick_width=1.5)
img.axis_type(1,"major",tick_len=6, tick_width=1.5)

dys = (-0.1, 0.3)
legend_pos = (0.02, 0.96)

col = 0

gauss_size = [2, 4, 2, 4]
sig_scale = [1.5, 1.5, 4, 4]
color = ["C0", "C4", "C1", "C3"]

for tag, filter_name in enumerate(sex_filters):
    print(filter_name, files)
    if tag < 2:
        for ii, sub_path in enumerate(data_path):
            h5f = h5py.File(data_path[ii] + "%s/%s/total.hdf5" % (filter_name, files), "r")
            data = h5f["/mc1"][()]
            h5f.close()
            # data = numpy.load(data_path + "%s/%s/total.npz" % (filter_name, files))
            lb = "$\\rm{MAG\_AUTO}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
            mc = data[:, ch]
            img.axs[0][ii].errorbar(x_coord, 100 * mc[0], 100 * mc[1], c=color[tag], linewidth=img.plt_line_width,
                        capsize=img.cap_size, label=lb, marker="s", fillstyle="none")

            h5f = h5py.File(data_path[ii] + "%s/flux2_ex1/total.hdf5" %filter_name, "r")
            data = h5f["/mc1"][()]
            h5f.close()
            # data = numpy.load(data_path + "%s/flux2_ex1/total.npz" % filter_name)
            lb = "$P_{k0}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
            mc = data[:, ch]
            img.axs[0][ii].errorbar(x_coord, 100 * mc[0], 100 * mc[1], c=color[tag + 2], linewidth=img.plt_line_width,
                        capsize=img.cap_size, label=lb, marker="p")
    else:
        for ii, sub_path in enumerate(data_path):
            h5f = h5py.File(data_path[ii] + "%s/%s/total.hdf5" % (filter_name, files), "r")
            data = h5f["/mc1"][()]
            h5f.close()
            # data = numpy.load(data_path + "%s/%s/total.npz" % (filter_name, files))
            lb = "$\\rm{MAG\_AUTO}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
            mc = data[:, ch]
            img.axs[0][ii].errorbar(x_coord, 100 * mc[0], 100 * mc[1], c=color[tag - 2], linewidth=img.plt_line_width,
                        capsize=img.cap_size, label=lb, marker="s", linestyle="--", alpha=0.65, fillstyle="none")


            h5f = h5py.File(data_path[ii] + "%s/flux2_ex1/total.hdf5" % filter_name, "r")
            data = h5f["/mc1"][()]
            h5f.close()
            # data = numpy.load(data_path + "%s/flux2_ex1/total.npz" % filter_name)
            lb = "$P_{k0}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
            mc = data[:, ch]
            img.axs[0][ii].errorbar(x_coord, 100 * mc[0], 100 * mc[1], c=color[tag], linewidth=img.plt_line_width,
                        capsize=img.cap_size, label=lb, marker="p", linestyle="--", alpha=0.65)

for i in range(2):
    xs = img.axs[0][i].set_xlim()
    ys = img.axs[0][i].set_ylim()
    img.axs[0][i].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width, c="grey", linestyle="--")
    img.axs[0][i].set_xlim(xs[0], xs[1])
    img.axs[0][i].set_ylim(ys[0] + dys[0], ys[1] + dys[1])
    img.axs[0][i].xaxis.set_major_formatter(xticks)
    img.axs[0][i].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
    img.axs[0][i].set_ylabel(ylabels[col], fontsize=img.xy_lb_size)
    img.axs[0][i].set_xticks(x_tick)
img.axs_text(0,0, 0.1, 0.65, "GII-sample",text_color="k", text_fontsize=img.xy_lb_size)
img.axs_text(0,1, 0.1, 0.65, "PII-sample",text_color="k", text_fontsize=img.xy_lb_size)
img.axs[0][0].legend(ncol=4, loc="upper left", fontsize=img.legend_size, bbox_to_anchor=(0.08,1.18))
img.save_img(pic_nm)
