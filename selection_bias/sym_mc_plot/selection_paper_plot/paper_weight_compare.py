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


data_path = "D:/cut/large_g/pts_bright/result"

ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]
ylabels = ["$c_1 \\times 10^4$", "$c_2 \\times 10^4$", "$c_1 \\times 10^4$", "$c_2 \\times 10^4$"]


x_tick = [i * cuts_num for i in range(0, ch_num,2)]
fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)


h5f = h5py.File(data_path + "/cuts_pi_all_sample_w_tflux_sq/sym/sex2_1.5/flux2_ex3/total.hdf5","r")
mc1_tf_pk_fit = h5f["/mc1"][()][:, ch]
mc2_tf_pk_fit = h5f["/mc2"][()][:, ch]
h5f.close()

h5f = h5py.File(data_path + "/cuts_pi_all_sample_w_tflux_sq/sym/sex2_1.5/mag_true/total.hdf5","r")
mc1_tf_tmag = h5f["/mc1"][()][:, ch]
mc2_tf_tmag = h5f["/mc2"][()][:, ch]
h5f.close()

h5f = h5py.File(data_path + "/cuts_pi_all_sample_w_maxpkfit_sq/sym/sex2_1.5/flux2_ex3/total.hdf5","r")
mc1_mpk_pk_fit = h5f["/mc1"][()][:, ch]
mc2_mpk_pk_fit = h5f["/mc2"][()][:, ch]
h5f.close()

h5f = h5py.File(data_path + "/cuts_pi_all_sample_w_maxpkfit_sq/sym/sex2_1.5/mag_true/total.hdf5","r")
mc1_mpk_tmag = h5f["/mc1"][()][:, ch]
mc2_mpk_tmag = h5f["/mc2"][()][:, ch]
h5f.close()


matplotlib.rcParams["font.family"] = "serif"
img = Image_Plot(fig_x=6,fig_y=4,xpad=0.2, ypad=0.2,legend_size=13)
img.subplots(1,1)
img.axis_type(0,"major",tick_len=6, tick_width=1.5)
img.axis_type(1,"major",tick_len=6, tick_width=1.5)

img.axs[0][0].errorbar(x_coord, 100 * mc1_tf_pk_fit[0], 100 * mc1_tf_pk_fit[1], linewidth=img.plt_line_width-0.5,
                       capsize=img.cap_size, marker="o", fillstyle="none",c="C2",label="$m_1$, weight=$F^{-2}$")
img.axs[0][0].errorbar(x_coord, 100 * mc2_tf_pk_fit[0], 100 * mc2_tf_pk_fit[1], linewidth=img.plt_line_width-0.5,
                       capsize=img.cap_size, marker="o", fillstyle="none", ls="--",c="C2",label="$m_2$, weight=$F^{-2}$")
# img.axs[0][0].errorbar(x_coord, 100 * mc1_tf_tmag[0], 100 * mc1_tf_tmag[1], linewidth=img.plt_line_width,
#                        capsize=img.cap_size, marker="s", fillstyle="none", ls="-.")
img.axs[0][0].errorbar(x_coord, 100 * mc1_mpk_pk_fit[0], 100 * mc1_mpk_pk_fit[1], linewidth=img.plt_line_width-0.5,
                       capsize=img.cap_size, marker="o", fillstyle="none",c="C1",label="$m_1$, weight=$\\tilde{\\nu}_F^{-2}$")
img.axs[0][0].errorbar(x_coord, 100 * mc2_mpk_pk_fit[0], 100 * mc2_mpk_pk_fit[1], linewidth=img.plt_line_width-0.5,
                       capsize=img.cap_size, marker="o", fillstyle="none", ls="--",c="C1",label="$m_2$, weight=$\\tilde{\\nu}_F^{-2}$")
# img.axs[0][0].errorbar(x_coord, 100 * mc1_mpk_tmag[0], 100 * mc1_mpk_tmag[1], linewidth=img.plt_line_width,
#                        capsize=img.cap_size, marker="s", fillstyle="none")
xs = img.axs[0][0].set_xlim()
ys = img.axs[0][0].set_ylim(-1.02, 1.02)
img.axs[0][0].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width-0.5, c="grey", alpha=0.7,linestyle="--")
img.axs[0][0].xaxis.set_major_formatter(xticks)
img.axs[0][0].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
img.axs[0][0].set_ylabel("$m_{1/2}\\times 10^2$", fontsize=img.xy_lb_size)
img.axs[0][0].set_xticks(x_tick)
img.axs[0][0].legend(ncol=2,fontsize=img.legend_size,frameon=False)
img.save_img("E:/weight_compare.pdf")
img.show_img()