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

source_b = "galsim_bright"
source_f = "galsim_dimmer"
# final pic name
pic_nm = "gal_mc.pdf"

# the bright source
total_path_1 = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/"%source_b
data_path_1 = total_path_1 + "cuts/sym/sex2_1.5/"
# the faint source
total_path_2 = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/"%source_f
data_path_2 = total_path_2 + "cuts/sym/sex2_1.5/"

names = ["P$_{k0}$", "MAG_AUTO", "MAG$_{true}$", "SNR$_S$", "SNR$_A$"]
files = ["flux2_ex1", "mag_auto", "flux2_ex5", "sex_snr", "snr_auto"]
colors = ["C0", "C2", "C1", "C3", "C4"]
ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]
ylabels = ["m$_1 \\times 10^2$", "m$_2 \\times 10^2$", "m$_1 \\times 10^2$", "m$_2 \\times 10^2$"]

fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

img = Image_Plot(fig_x=8, fig_y=5)
img.subplots(2,2)
# img.set_style_default()
# img.set_style()
# dys = [(-0.2, 1), (-0.2, 1), (-0.4, 0.5), (-0.4, 0.5)]
dys = [(-0.1, 0.35), (-0.1, 0.3), (0, 0.1), (0, 0.1)]

for j in range(4):

    row, col = divmod(j, 2)
    if row == 0:
        h5f_all = h5py.File(total_path_1 + "data/shear_result.hdf5", "r")
    else:
        h5f_all = h5py.File(total_path_2 + "data/shear_result.hdf5", "r")
    mc_all = h5f_all["/mc%d" % (col + 1)].value
    h5f_all.close()

    img.axs[row][col].errorbar(x_coord[0], mc_all[0]*100, mc_all[1]*100, c="grey", capsize=img.cap_size,
                               marker="o", mfc="none", label="All sources",linewidth=img.plt_line_width)
    for i in range(len(files)):
        if row == 0:
            #             bright source
            # data = numpy.load(data_path_1 + files[i] + "/total.npz")
            h5f = h5py.File(data_path_1 + files[i] + "/total.hdf5","r")
        else:
            #             faint source
            # data = numpy.load(data_path_2 + files[i] + "/total.npz")
            h5f = h5py.File(data_path_2 + files[i] + "/total.hdf5", "r")

        mc = h5f["/mc%d"%(col+1)].value[:,ch]
        h5f.close()

        img.axs[row][col].errorbar(x_coord, 100 * mc[0], 100 * mc[1], c=colors[i], linewidth=img.plt_line_width,
                    capsize=img.cap_size, label=names[i], marker="o",mfc="none")

    xs = img.axs[row][col].set_xlim()
    ys = img.axs[row][col].set_ylim()
    img.axs[row][col].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width, c="grey", linestyle="--")
    img.axs[row][col].set_xlim(xs[0], xs[1])

    img.axs[row][col].set_ylim(ys[0] + dys[j][0], ys[1] + dys[j][1])

    img.axs[row][col].xaxis.set_major_formatter(xticks)
    img.axs[row][col].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
    img.axs[row][col].set_ylabel(ylabels[j], fontsize=img.xy_lb_size)
    img.axs[row][col].legend(ncol=2, loc="best", fontsize=img.legend_size, frameon=False)#, bbox_to_anchor=legend_pos[j]
img.figure.subplots_adjust(hspace=0.25, wspace=0.25)
img.save_img(pic_nm)


