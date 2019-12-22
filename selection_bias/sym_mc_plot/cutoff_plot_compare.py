import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
path.append("D:/Github/astrophy-research/mylib")
from sys import argv
import h5py
from plot_tool import Image_Plot
import matplotlib.ticker as mtick


# the number of thread must be equal to cutoff number
total_path = argv[1]
filter_name = argv[2]

files = ["mag_true", "mag_auto", "rfactor", "flux2_ex1", "flux2_ex2", "flux2_ex3", "snr_sex"]
file_legend = ["mag_true", "mag_auto", "resolution", "Pk0_1", "Pk0_2", "Pk0_3", "snr_sex"]
colors = ["C%d"%i for i in range(len(files))]

fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)
ylabels = ["$m_1[10^{-2}]$", "$c_1[10^{-4}]$","$m_2[10^{-2}]$", "$c_2[10^{-4}]$"]

pic_path = total_path + "/%s/compare_%s.png"%(filter_name, filter_name)
ch_num = 10
x_coord = [i * 10 for i in range(ch_num)]
ch = [i for i in range(ch_num)]
m_scale, c_scale = 100, 10000
img = Image_Plot(fig_x=6,fig_y=4, xpad=0.25, ypad=0.25,plt_line_width=1.5)
img.subplots(2,2)

for tag, nm in enumerate(files):
    file_path = total_path + "/%s/%s/total.hdf5"%(filter_name,nm)
    try:
        h5f = h5py.File(file_path,"r")
        mc1 = h5f["/mc1"][()]
        mc2 = h5f["/mc2"][()]
        h5f.close()

        img.axs[0][0].errorbar(x_coord, mc1[0]*m_scale, mc1[1]*m_scale, c=colors[tag], linewidth=img.plt_line_width, capsize=img.cap_size,
                               label=file_legend[tag], marker="s")
        img.axs[1][0].errorbar(x_coord, mc2[0]*m_scale, mc2[1]*m_scale, c=colors[tag], linewidth=img.plt_line_width, capsize=img.cap_size,
                               label=file_legend[tag], marker="s")

        img.axs[0][1].errorbar(x_coord, mc1[2]*c_scale, mc1[3]*c_scale, c=colors[tag], linewidth=img.plt_line_width, capsize=img.cap_size,
                               label=file_legend[tag], marker="s")
        img.axs[1][1].errorbar(x_coord, mc2[2]*c_scale, mc2[3]*c_scale, c=colors[tag], linewidth=img.plt_line_width, capsize=img.cap_size,
                               label=file_legend[tag], marker="s")

    except:
        print("Can't find the file!", file_path)

for i in range(2):
    for j in range(2):
        xs = img.axs[i][j].set_xlim()
        img.axs[i][j].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width, c="grey", linestyle="--")
        img.axs[i][j].set_xlim(xs[0], xs[1])
        img.axs[i][j].xaxis.set_major_formatter(xticks)
        img.axs[i][j].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
        tag = 2*i+j
        img.axs[i][j].set_ylabel(ylabels[tag], fontsize=img.xy_lb_size)

img.axs[0][0].legend(fontsize=img.legend_size-2,ncol=len(files),loc="lower left", bbox_to_anchor=(0.01,1.02))
img.save_img(pic_path)
img.show_img()