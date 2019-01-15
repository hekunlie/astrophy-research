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
xticks = mtick.FormatStrFormatter(fmt)
fig = plt.figure(figsize=figs)
ax = fig.add_subplot(111)
ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
for axis in ["bottom", "left", "top", "right"]:
    # the line width of the frame
    ax.spines[axis].set_linewidth(axis_linewidth)
ax.xaxis.set_tick_params(which="both",direction="in",length=6, width=axis_linewidth)

total_path = "/mnt/ddnfs/data_users/hkli/simu_test/"

para_path = total_path + "parameters/para_0.hdf5"
para_h5 = h5py.File(para_path, "r")
mag_t = para_h5["/mag"].value
para_h5.close()

fourier_path = total_path +"result/data/data_1.5sig/data_0.hdf5"
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
flux_err[detect_s] = 1
snr_auto = flux_auto / flux_err
scale_sa = numpy.load(total_path+"result/cuts/cache_snr_auto.npz")["arr_1"][0]
# snr
snr = s_data[:, 0]
scale_s = numpy.load(total_path+"result/cuts/cache_sex_snr.npz")["arr_1"][0]
# mag_auto
mag = s_data[:, 3]
scale_m = numpy.load(total_path+"result/cuts/cache_mag_auto.npz")["arr_1"][0]

num, bins = ax.hist(mag_t, 60, label="TOTAL SAMPLE")[:2]
ax.hist(mag_t[detect_f], bins, label="Fourier Quad")
ax.hist(mag_t[detect_s], bins, label="SExtractor")



