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



matplotlib.rcParams["font.family"] = "serif"

#  the m & c of different filters from point source
source_b = "pts_ellip_psf"
pic_nm = "pts_filters_f_c2.pdf"

# the bright source
data_path_1 = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/cuts/sym/"  % source_b

sex_filters = ["sex2_1.5", "sex4_1.5", "sex2_4", "sex4_4"]
names = ["P$_{k0}$", "MAG_AUTO", "MAG$_{true}$", "SNR$_S$", "SNR$_A$"]
files = "mag_auto"
ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]
ylabels = ["$c_1 \\times 10^4$", "$c_2 \\times 10^4$", "$c_1 \\times 10^4$", "$c_2 \\times 10^4$"]


x_tick = [i * cuts_num for i in range(0,ch_num,2)]
fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

img = Image_Plot(fig_x=9,fig_y=6)
img.subplots(1,1)

dys = (-0.1, 0.46)
legend_pos = (0.02, 0.01)
plt_line_width = 2
col = 0


gauss_size = [2, 4, 2, 4]
sig_scale = [1.5, 1.5, 4, 4]
color = ["C0", "C4", "C1", "C3"]

for tag, filter_name in enumerate(sex_filters):
    print(filter_name, files)
    if tag < 2:
        h5f = h5py.File(data_path_1 + "%s/%s/total.hdf5" % (filter_name, files))
        lb = "$MAG\_AUTO, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = h5f['/mc2'].value[:, ch]
        img.axs[0][0].errorbar(x_coord, 10000 * mc[2], 10000 * mc[3], c=color[tag], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=lb, marker="s", fillstyle="none")
        h5f.close()

        h5f = h5py.File(data_path_1 + "%s/flux2_ex1/total.hdf5" % filter_name)
        lb = "$P_{k0}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = h5f['/mc2'].value[:, ch]
        img.axs[0][0].errorbar(x_coord, 10000 * mc[2], 10000 * mc[3], c=color[tag + 2], linewidth=img.plt_line_width,
                               capsize=img.cap_size, label=lb, marker="p")
        h5f.close()
    else:
        h5f = h5py.File(data_path_1 + "%s/%s/total.hdf5" % (filter_name, files))
        lb = "$MAG\_AUTO, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = h5f['/mc2'].value[:, ch]
        img.axs[0][0].errorbar(x_coord, 10000 * mc[2], 10000 * mc[3], c=color[tag - 2], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=lb, marker="s", linestyle="--", alpha=0.6, fillstyle="none")
        h5f.close()

        h5f = h5py.File(data_path_1 + "%s/flux2_ex1/total.hdf5" % filter_name)
        lb = "$P_{k0}, gauss\_{%.1f}, %.1f\sigma$" % (gauss_size[tag], sig_scale[tag])
        mc = h5f['/mc2'].value[:, ch]
        img.axs[0][0].errorbar(x_coord, 10000 * mc[2], 10000 * mc[3], c=color[tag], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=lb, marker="p", linestyle="--", alpha=0.6)
        h5f.close()
xs = img.axs[0][0].set_xlim()
ys = img.axs[0][0].set_ylim(-3.1, 4)
img.axs[0][0].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width, c="grey", linestyle="--")
img.axs[0][0].set_xlim(xs[0], xs[1])
# ax.set_ylim(ys[0] + dys[0], ys[1]+ dys[1])
img.axs[0][0].xaxis.set_major_formatter(xticks)
img.axs[0][0].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
img.axs[0][0].set_ylabel(ylabels[3], fontsize=img.xy_lb_size)
img.axs[0][0].set_xticks(x_tick)
img.axs[0][0].legend(ncol=2, loc="lower left", fontsize=img.legend_size-2, bbox_to_anchor=legend_pos, frameon=False)
img.save_img(pic_nm)
