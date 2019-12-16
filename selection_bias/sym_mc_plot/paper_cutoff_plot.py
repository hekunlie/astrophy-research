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
pic_nm = "gal_mc_original_pdf.pdf"
pic_nm_png = pic_nm.split(".")[0]+".png"
# the file name
file_name = "cuts_bk_2019_12_9_original_pdf/"
sex_filter_name = "sym/sex2_1.5/"
# result file of all source
shear_result_all = file_name + "shear_result.hdf5"
# the bright source
total_path_1 = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/"%source_b
data_path_1 = total_path_1 + file_name + sex_filter_name
# the faint source
total_path_2 = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/"%source_f
data_path_2 = total_path_2 + file_name + sex_filter_name

names = ["P$_{k0}$", "MAG_AUTO", "MAG$_{true}$", "SNR$_S$", "Resolution"]
files = ["flux2_ex1", "mag_auto", "mag_true", "snr_sex", "rfactor"]
colors = ["C0", "C2", "C1", "C3", "C4"]
ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]
ylabels = ["m$_1 \\times 10^2$", "m$_2 \\times 10^2$", "m$_1 \\times 10^2$", "m$_2 \\times 10^2$"]

fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

img = Image_Plot(fig_x=6, fig_y=4, ypad=0.2, xpad=0.2)
img.subplots(2,2)
img.axis_type(0,"major",tick_len=8, tick_width=2)
img.axis_type(1,"major",tick_len=8, tick_width=2)


# # pts sample
# text_pos = [[6, 1.],[6, 1.],[6, 1.6],[6, 1.6]]
# sample_name = ["PI sample","PI sample","PII sample","PII sample"]
# xy_lims = [(-2.9, 1.4),(-2.9, 1.8),(-4.98, 2.6),(-4.9, 2.6)]

# galsim sample
text_pos = [[6, 0.7],[6, 0.7],[6, 0.7],[6, 0.7]]
sample_name = ["GI sample","GI sample","GII sample","GII sample"]
xy_lims = [(-0.85, 0.9),(-0.91, 0.9),(-2, 1.1),(-2, 1.1)]

for j in range(4):

    row, col = divmod(j, 2)
    if row == 0:
        h5f_all = h5py.File(total_path_1 + shear_result_all, "r")
    else:
        h5f_all = h5py.File(total_path_2 + shear_result_all, "r")
    mc_all = h5f_all["/mc%d" % (col + 1)][()]
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

        mc = h5f["/mc%d"%(col+1)][()][:,ch]
        h5f.close()

        img.axs[row][col].errorbar(x_coord, 100 * mc[0], 100 * mc[1], c=colors[i], linewidth=img.plt_line_width-0.25,
                    capsize=img.cap_size, label=names[i], marker="d",mfc="none")

    xs = img.axs[row][col].set_xlim()
    ys = img.axs[row][col].set_ylim()
    img.axs[row][col].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width, c="grey", linestyle="--")
    img.axs[row][col].set_xlim(xs[0], xs[1])

    # img.axs[row][col].set_ylim(ys[0] + dys[j][0], ys[1] + dys[j][1])
    # img.axs[row][col].set_ylim((-6.5,4.5))

    img.axs[row][col].xaxis.set_major_formatter(xticks)
    img.axs[row][col].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
    img.axs[row][col].set_ylabel(ylabels[j], fontsize=img.xy_lb_size)
    # img.axs[row][col].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
    # if j == 0 or j == 2:
    #     img.axs[row][col].set_ylabel("m$\\times 10^2$", fontsize=img.xy_lb_size)
    #     img.axs[row][col].set_yticklabels([])

    if j == 0:
        img.axs[row][col].legend(ncol=6, fontsize=img.legend_size-2,loc='upper left', bbox_to_anchor=(0.11,1.17))

    img.axs[row][col].set_ylim(xy_lims[j])
    img.axs[row][col].text(text_pos[j][0], text_pos[j][1], sample_name[j], color='k', ha='left', va='center',
                           fontsize=img.legend_size - 1,fontweight="medium")

    # if j ==2 or j == 3:
    #     img.axs[row][col].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
    # else:
    #     img.axs[row][col].set_xticklabels([])


# img.figure.subplots_adjust(hspace=0.25, wspace=0.25)
img.save_img(pic_nm)
img.save_img(pic_nm_png)


