import matplotlib
# matplotlib.use("Agg")
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

mc_cmd = 0

source_b = "pts_bright"
source_f = "pts_dimmer"
parent_path = "D:/cut/large_g"
# final pic name
pic_nm = "/pts_mc.pdf"
pic_nm_png = pic_nm.split(".")[0]+".png"
# the file name
file_name = "cuts_pi_half_sample/"
sex_filter_name = "sym/sex2_1.5/"
# result file of all source
shear_result_all = file_name + "shear_result.hdf5"
# the bright source
total_path_1 = parent_path + "/%s/result/"%source_b
data_path_1 = total_path_1 + file_name + sex_filter_name
# the faint source
total_path_2 = parent_path + "/%s/result/"%source_f
data_path_2 = total_path_2 + file_name + sex_filter_name

names = ["$\\nu_F$", "MAG_AUTO", "MAG$_{\\rm{true}}$", "SNR", "Resolution"]
files = ["flux2_ex1", "mag_auto", "mag_true", "snr_sex", "rfactor"]
colors = ["C0", "C2", "C1", "C3", "C4"]
ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]


fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

img = Image_Plot(fig_x=6, fig_y=4, ypad=0.25, xpad=0.25,pts_size=5,cap_size=3)
img.subplots(2,2)
img.axis_type(0,"major",tick_len=8, tick_width=1.5)
img.axis_type(1,"major",tick_len=8, tick_width=1.5)


# # pts sample
if mc_cmd == 0:
    ylabels = ["m$_1 \\times 10^2$", "m$_2 \\times 10^2$", "m$_1 \\times 10^2$", "m$_2 \\times 10^2$"]
    text_pos = [[50, 2.5],[50, 2.5],[50, 2.2],[50, 2.2]]
    sample_name = ["PI sample","PI sample","PII sample","PII sample"]
    xy_lims = [(-4.1, 3.7),(-4.1, 3.7),(-6.5, 3.7),(-6.5, 3.7)]


# # galsim sample
elif mc_cmd == 1:
    ylabels = ["m$_1 \\times 10^2$", "m$_2 \\times 10^2$", "m$_1 \\times 10^2$", "m$_2 \\times 10^2$"]
    text_pos = [[6, 1.6],[6, 1.6],[6, 1.4],[6, 1.4]]
    sample_name = ["GI sample","GI sample","GII sample","GII sample"]
    xy_lims = [(-1.9, 2.2),(-1.9, 2.2),(-3.6, 2.2),(-3.6, 2.2)]

# epsf sample
else:
    ylabels = ["c$_1 \\times 10^3$", "c$_2 \\times 10^3$", "c$_1 \\times 10^3$", "c$_2 \\times 10^3$"]
    text_pos = [[3, 0.075],[55, 0.25],[3, 0.065],[55, 0.6]]
    sample_name = ["GI sample","GI sample","PI sample","PI sample"]
    xy_lims = [(-0.042, 0.092),(-0.33, 0.33),(-0.12,0.092),(-0.71, 0.78)]

for j in range(4):

    row, col = divmod(j, 2)
    if row == 0:
        h5f_all = h5py.File(total_path_1 + shear_result_all, "r")
    else:
        h5f_all = h5py.File(total_path_2 + shear_result_all, "r")
    mc_all = h5f_all["/mc%d" % (col + 1)][()]
    h5f_all.close()
    if mc_cmd < 2:
        img.axs[row][col].errorbar(x_coord[0]-1.5, mc_all[0]*100, mc_all[1]*100, c="dimgray", capsize=img.cap_size,
                                   marker="o", label="All sources",linewidth=img.plt_line_width-0.5, ms=img.pts_size)
    else:
        img.axs[row][col].errorbar(x_coord[0]-1.5, mc_all[2]*1000, mc_all[3]*1000, c="dimgrey", capsize=img.cap_size,
                                   marker="o", label="All sources",linewidth=img.plt_line_width-0.5,ms=img.pts_size)
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
        if mc_cmd < 2:
            img.axs[row][col].errorbar(x_coord, 100 * mc[0], 100 * mc[1], c=colors[i], linewidth=img.plt_line_width-0.5,
                        capsize=img.cap_size, label=names[i], marker="o",mfc="none",ms=img.pts_size)
        else:
            img.axs[row][col].errorbar(x_coord, 1000 * mc[2], 1000 * mc[3], c=colors[i], linewidth=img.plt_line_width-0.5,
                        capsize=img.cap_size, label=names[i], marker="o",mfc="none",ms=img.pts_size)

    xs = img.axs[row][col].set_xlim()
    ys = img.axs[row][col].set_ylim()
    img.axs[row][col].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width-0.75, c="grey", linestyle="--",alpha=0.75)
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
        img.axs[row][col].legend(ncol=6, fontsize=img.legend_size-2,
                                 loc='upper left', bbox_to_anchor=(0.11,1.18), edgecolor="k",fancybox=False)

    img.axs[row][col].set_ylim(xy_lims[j])
    img.axs[row][col].text(text_pos[j][0], text_pos[j][1], sample_name[j], color='k', ha='left', va='center',
                           fontsize=img.legend_size - 1,fontweight="medium")

    # if j ==2 or j == 3:
    #     img.axs[row][col].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
    # else:
    #     img.axs[row][col].set_xticklabels([])


# img.figure.subplots_adjust(hspace=0.25, wspace=0.25)

# img.save_img(parent_path + pic_nm)
# img.save_img(parent_path + pic_nm_png)
img.save_img("E:" + pic_nm)
img.save_img("E:" + pic_nm_png)

img.show_img()


exit()
matplotlib.rcParams["font.family"] = "serif"

mc_cmd = 2

source_b = "galsim_bright_epsf"
source_f = "pts_bright_epsf"

parent_path = "D:/cut/large_g"
# final pic name
pic_nm = "/epsf_1_mc.pdf"
pic_nm_png = pic_nm.split(".")[0]+".png"
# the file name
file_name = "cuts_pi_half_sample/"
sex_filter_name = "sym/sex2_1.5/"
# result file of all source
shear_result_all = file_name + "shear_result.hdf5"
# the bright source
total_path_1 = parent_path + "/%s/result/"%source_b
data_path_1 = total_path_1 + file_name + sex_filter_name
# the faint source
total_path_2 = parent_path + "/%s/result/"%source_f
data_path_2 = total_path_2 + file_name + sex_filter_name

# the bright source
total_path_1_ = parent_path + "/%s_1/result/"%source_b
data_path_1_ = total_path_1_ + file_name + sex_filter_name
# the faint source
total_path_2_ = parent_path + "/%s_1/result/"%source_f
data_path_2_ = total_path_2_ + file_name + sex_filter_name

names = ["$\\nu_F$", "MAG_AUTO", "MAG$_{\\rm{true}}$", "SNR", "Resolution"]
files = ["flux2_ex1", "mag_auto", "mag_true", "snr_sex", "rfactor"]
colors = ["C0", "C2", "C1", "C3", "C4"]
ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]


fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

img = Image_Plot(fig_x=5, fig_y=3, ypad=0.3, xpad=0.3,pts_size=5,cap_size=3)
img.subplots(2,3)
img.axis_type(0,"major",tick_len=8, tick_width=1.5)
img.axis_type(1,"major",tick_len=8, tick_width=1.5)


ylabels = ["c$_1 \\times 10^3$", "c$_2 \\times 10^3$",
           "c$_2 \\times 10^3$", "c$_1 \\times 10^3$",
           "c$_2 \\times 10^3$", "c$_2 \\times 10^3$"]

text_pos = [[3, 0.075],[55, 0.25],[55, 0.25],
            [3, 0.067],[55, 0.6],[55, 0.6]]

sample_name = ["GI sample","GI sample","GI sample",
               "PI sample", "PI sample","PI sample"]
xy_lims = [(-0.042, 0.092),(-0.33, 0.33),(-0.33, 0.33),
           (-0.12,0.092),(-0.71, 0.78),(-0.71, 0.78)]

c12_label = [[1,2,2],[1,2,2]]
for j in range(6):

    row, col = divmod(j, 3)
    if row == 0:
        if col <2:
            h5f_all = h5py.File(total_path_1 + shear_result_all, "r")
            print(row, col, 1)
        else:
            h5f_all = h5py.File(total_path_1_ + shear_result_all, "r")
            print(row, col, 2)
    else:
        if col<2:
            h5f_all = h5py.File(total_path_2 + shear_result_all, "r")
            print(row, col, 1)
        else:
            h5f_all = h5py.File(total_path_2_ + shear_result_all, "r")
            print(row, col, 2)
    mc_all = h5f_all["/mc%d" %c12_label[row][col]][()]
    h5f_all.close()

    img.axs[row][col].errorbar(x_coord[0]-1.5, mc_all[2]*1000, mc_all[3]*1000, c="dimgrey", capsize=img.cap_size,
                                   marker="o", label="All sources",linewidth=img.plt_line_width-0.5,ms=img.pts_size)
    for i in range(len(files)):
        if row == 0:
            if col< 2:
                #             bright source
                # data = numpy.load(data_path_1 + files[i] + "/total.npz")
                h5f = h5py.File(data_path_1 + files[i] + "/total.hdf5","r")
            else:
                h5f = h5py.File(data_path_1_ + files[i] + "/total.hdf5","r")

        else:
            if col<2:
                #             faint source
                # data = numpy.load(data_path_2 + files[i] + "/total.npz")
                h5f = h5py.File(data_path_2 + files[i] + "/total.hdf5", "r")
            else:
                h5f = h5py.File(data_path_2_ + files[i] + "/total.hdf5", "r")

        mc = h5f["/mc%d"%c12_label[row][col]][()][:,ch]
        h5f.close()

        img.axs[row][col].errorbar(x_coord, 1000 * mc[2], 1000 * mc[3], c=colors[i], linewidth=img.plt_line_width-0.5,
                        capsize=img.cap_size, label=names[i], marker="o",mfc="none",ms=img.pts_size)

    xs = img.axs[row][col].set_xlim()
    ys = img.axs[row][col].set_ylim()
    img.axs[row][col].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width-0.75, c="grey", linestyle="--",alpha=0.75)
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
        img.axs[row][col].legend(ncol=6, fontsize=img.legend_size-2,
                                 loc='upper left', bbox_to_anchor=(0.57,1.23), edgecolor="k",fancybox=False)

    img.axs[row][col].set_ylim(xy_lims[j])
    img.axs[row][col].text(text_pos[j][0], text_pos[j][1], sample_name[j], color='k', ha='left', va='center',
                           fontsize=img.legend_size - 1,fontweight="medium")

    # if j ==2 or j == 3:
    #     img.axs[row][col].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
    # else:
    #     img.axs[row][col].set_xticklabels([])


# img.figure.subplots_adjust(hspace=0.25, wspace=0.25)

# img.save_img(parent_path + pic_nm)
# img.save_img(parent_path + pic_nm_png)
img.save_img("E:" + pic_nm)
img.save_img("E:" + pic_nm_png)

img.show_img()