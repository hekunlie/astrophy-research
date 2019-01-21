import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import ScalarFormatter
import tool_box
import h5py


fmt='%2.f%%'
fig_x = 8
fig_y = fig_x*4/6
figs = (fig_x, fig_y)
fonts = 20
xy_lb_size = 18
legend_size = fonts - 4
axis_linewidth = 1.2
plt_line_width = 2
cap_size = 5
# xticks = mtick.FormatStrFormatter(fmt)
fig = plt.figure(figsize=figs)
ax = fig.add_subplot(111)
ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
for axis in ["bottom", "left", "top", "right"]:
    # the line width of the frame
    ax.spines[axis].set_linewidth(axis_linewidth)
ax.xaxis.set_tick_params(which="both",direction="in",length=6, width=axis_linewidth)

total_path = "/mnt/ddnfs/data_users/hkli/simu_test1/"

para_path = total_path + "parameters/para_0.hdf5"
para_h5 = h5py.File(para_path, "r")
mag_t = para_h5["/mag"].value
para_h5.close()

fourier_path = total_path + "result/data/data_1.5sig/data_0.hdf5"
f_h5 = h5py.File(fourier_path, "r")
f_data = f_h5["/data"].value
f_h5.close()

sex_path = total_path + "result/data/sex2_1.5/sex_0.npz"
s_data = numpy.load(sex_path)["arr_0"]

detect_f = f_data[:, -1] > -1
detect_s = s_data[:, 0] > 0

# P(k=0)
pk0 = f_data[:, 4]
scale_p = numpy.load(total_path+"result/cuts/cache_flux2_ex1.npz")["arr_1"][0]

# snr_auto
flux_auto = s_data[:, 1]
flux_err = s_data[:, 2]
idx_s = flux_err == 0
flux_err[idx_s] = 1
snr_auto = flux_auto / flux_err
scale_sa = numpy.load(total_path+"result/cuts/cache_snr_auto.npz")["arr_1"][0]
# snr
snr = s_data[:, 0]
scale_s = numpy.load(total_path+"result/cuts/cache_sex_snr.npz")["arr_1"][0]
# mag_auto
mag = s_data[:, 3]
scale_m = -numpy.load(total_path+"result/cuts/cache_mag_auto.npz")["arr_1"][0]

cut_data = [pk0, snr_auto, snr, mag]
cut_data_scales = [scale_p, scale_sa, scale_s,  scale_m]

num, bins = ax.hist(mag_t, 60, label="TOTAL SAMPLE")[:2]
ax.hist(mag_t[detect_f], bins, label="Fourier Quad")
ax.hist(mag_t[detect_s], bins, label="SExtractor")
xs = ax.set_xlim()
ys = ax.set_ylim()

labels = ["P(k=0)", "SNR_AUTO", "SNR", "MAG_AUTO"]
cuts_num = 10
for i in range(len(cut_data)):
    if i == 0:
        idx = detect_f
    else:
        idx = detect_s
    if i == 3:
        cut_data_sort = numpy.sort(-cut_data[i][idx])
        cut_data_step = int(len(cut_data_sort) / cuts_num)
        cutoff_scalar = [-cut_data_sort[i * cut_data_step] for i in range(cuts_num)]
    else:
        cut_data_sort = numpy.sort(cut_data[i][idx])
        cut_data_step = int(len(cut_data_sort) / cuts_num)
        cutoff_scalar = [cut_data_sort[i * cut_data_step] for i in range(cuts_num)]
    print(labels[i], cutoff_scalar)
    for j in range(2):
        cut_s = cutoff_scalar[4+j*3]
        tag = numpy.where(cut_s == cut_data[i][idx])
        # mag_refer = mag_t[idx][tag]

        # print(labels[i], cut_s, tag)
        # if j == 0:
        #     lines = "-"
        #     lb = labels[i] + " 30%"
        # else:
        #     lines = "--"
        #     lb = labels[i] + " 60%"
        # ax.plot([mag_refer, mag_refer], [ys[0], ys[1]], c="C%d"%i, linestyle=lines, label=lb)
# ax.legend()
# ax.set_ylim(ys[0], ys[1])
# plt.savefig("/home/hkli/work/mag_cut.pdf",bbox_inches='tight')