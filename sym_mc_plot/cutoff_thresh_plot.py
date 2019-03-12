import numpy
import h5py
import matplotlib.pyplot as plt

file_no = 0
# data of Fourier Quad
pk0_scale = numpy.load("F:/cache_flux2_ex1.npz")["arr_1"][file_no]/64/60
h5f = h5py.File("F:/data_%d.hdf5"%file_no)
data = h5f["/data"].value
h5f.close()
idx = data[:, -1] > 0
pk0 = data[:, 4][idx]/64/60

# data of SExtractor
snr_idx = 0
flux_auto_idx = 1
flux_err_idx = 2
mag_auto_idx = 3

sex_data = numpy.load("F:/sex_%d.npz"%file_no)["arr_0"]
idx = sex_data[:,0] > 0

snr_scale = numpy.load("F:/cache_sex_snr.npz")["arr_1"][file_no]
snr = sex_data[:, snr_idx][idx]

mag_scale = numpy.load("F:/cache_mag_auto.npz")["arr_1"][file_no]
mag = sex_data[:,mag_auto_idx][idx]

snr_auto_scale = numpy.load("F:/cache_snr_auto.npz")["arr_1"][file_no]
flux_auto = sex_data[:, flux_auto_idx][idx]
flux_err = sex_data[:, flux_err_idx][idx]
snr_auto = flux_auto/flux_err

scales = [pk0_scale, snr_scale, mag_scale, snr_auto_scale]
datas = [pk0, snr, mag, snr_auto]
names = ["P$_{k0}$", "SNR$_S$", "MAG_AUTO", "SNR$_A$"]

fmt = '%2.f%%'
fig_x = 8
fig_y = fig_x * 4 / 6
figs = (fig_x, fig_y)
fonts = 20
xy_lb_size = 18
legend_size = fonts - 4
axis_linewidth = 1.2
plt_line_width = 2
cap_size = 5

i = 1

x_max = scales[i][8]+(scales[i][8] - scales[i][7])*0.5
idx = datas[i] < x_max
fig = plt.figure(figsize=figs)
ax = fig.add_subplot(111)
ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
for axis in ["bottom", "left", "top", "right"]:
    # the line width of the frame
    ax.spines[axis].set_linewidth(axis_linewidth)
ax.xaxis.set_tick_params(which="both", direction="in", length=5, width=axis_linewidth)
ax.yaxis.set_tick_params(which="major", direction="in", length=5, width=axis_linewidth)
ax.yaxis.set_tick_params(which="minor", direction="in", length=5, width=axis_linewidth)
#ax.yaxis.set_powerlimits((-1,1))
ax.set_xlabel(names[i],fontsize=fonts)
ax.set_ylabel("Galaxy number",fontsize=fonts)

ax.hist(datas[i][idx],50)

ys = ax.set_ylim()
xs = ax.set_xlim()
for ic in range(len(scales[i])-1):
    ax.plot([scales[i][ic], scales[i][ic]],[0, ys[1]], linestyle="--", c="C1", linewidth=2)
ax.set_ylim(ys[0], ys[1])
ax.set_xlim(xs[0], x_max)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.show()
plt.close()