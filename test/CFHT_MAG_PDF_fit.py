import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import ScalarFormatter
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import numpy
import tool_box

data_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/chft_cat/"
cat_path = data_path + "cat.npz"
if os.path.exists(cat_path):
    data = numpy.load(cat_path)["arr_0"]
else:
    files = os.listdir(data_path)
    t = 0
    for data_f in files:
        if".dat" in data_f:
            print(data_f)
            temp = numpy.loadtxt(data_path + data_f, skiprows=1)
            if t == 0:
                data = temp
            else:
                data = numpy.row_stack((data,temp))
            t+=1
    numpy.savez(cat_path, data)

cfht_mag = data[:,-1]
print(cfht_mag.shape, cfht_mag.max(), cfht_mag.min())

idx_1 = cfht_mag > 0
idx_2 = cfht_mag < 50
mag = cfht_mag[idx_1&idx_2]

bin_num = 50
nums, bins = numpy.histogram(mag, bin_num)
print(nums, nums.shape)
print(bins, bins.shape)

fit_x = bins[:bin_num]
print(fit_x, fit_x.shape)
idx_1 = fit_x <= 24
idx_2 = fit_x >= 18.5
fit_nums = numpy.log10(nums[idx_1&idx_2])
print(fit_nums.shape)
paras = tool_box.fit_1d(fit_x[idx_1&idx_2],fit_nums, 1, "lsq")
print(paras)

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
ax.yaxis.set_tick_params(which="major",direction="in",length=8, width=axis_linewidth)
ax.yaxis.set_tick_params(which="minor",direction="in",length=4, width=axis_linewidth)

x = numpy.linspace(14, 28, 500)
fxy = 10**(paras[0] + paras[1]*x)
ax.scatter(bins[:bin_num],nums, color="black", label="CFHTLenS")
ax.plot(x, fxy, color="C1", label="Fit")
ax.set_xlim(19, 27.2)
ys = ax.set_ylim()
ax.set_ylim(6000, ys[1])
#x1, x2 = ax.set_xlim(18.9,25.1)

ax.set_yscale("log")
ax.set_ylabel("Galaxy number",fontsize=xy_lb_size)
ax.set_xlabel("Magnitude",fontsize=xy_lb_size)
ax.legend(fontsize=xy_lb_size)
plt.savefig("mag_fit.pdf",bbox_inches='tight')
plt.show()