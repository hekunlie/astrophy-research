import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import ScalarFormatter
from sys import path
path.append('D:/Github/astrophy-research/my_lib')
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
data_path = "F:\works\selection_bias\point\\100_point_20-25\sym\sex2_1.5/"
data_path = data_path.replace("\\","/")
names = ["P$_{k0}$", "MAG_AUTO", "MAG$_{true}$", "SNR", "SNR_AUTO"]
files = ["flux2_ex1","mag_auto", "flux2_ex5","sex_snr","snr_auto"]

ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]

plt_item = 3
ylabels = ["m$_1 \\times 10^2$", "m$_2 \\times 10^2$",
           "c$_1 \\times 10^4$", "c$_2 \\times 10^4$"]
mc_item = ["m1", "m2", "c1", "c2"]
pic_name = mc_item[plt_item] + "_pts_b.pdf"
for i in range(len(files)):
    data = numpy.load(data_path+files[i]+"/total.npz")
    mc1 = data['arr_0'][:, ch]
    mc2 = data['arr_1'][:, ch]
    if i == 2:
        line_style ="-"
    else:
        line_style = "-"
    if plt_item == 0:
        # m1
        ax.errorbar(x_coord, 100*(mc1[0] - 1), 100*mc1[1], c="C%d"%i,linewidth=plt_line_width,
                    capsize=cap_size, label=names[i], marker="s")
    elif plt_item == 1:
        # m2
        ax.errorbar(x_coord, 100*(mc2[0] - 1), 100*mc2[1], c="C%d"%i,linewidth=plt_line_width,
                    capsize=cap_size, label=names[i], marker="s")
    elif plt_item == 2:
        # c1
        ax.errorbar(x_coord, 10000*mc1[2], 10000*mc1[3], c="C%d"%i,linewidth=plt_line_width,
                    capsize=cap_size, label=names[i], marker="s", linestyle=line_style)
    else:
        # c2
        ax.errorbar(x_coord, 10000*mc2[2], 10000*mc2[3], c="C%d"%i,linewidth=plt_line_width,
                    capsize=cap_size, label=names[i], marker="s")
xs = ax.set_xlim()
ys = ax.set_ylim()
ax.plot([xs[0],100],[0,0], linewidth=plt_line_width, c="grey", linestyle="--")
ax.set_xlim(xs[0], xs[1])
ax.set_ylim(ys[0], ys[1]+0.3)
ax.xaxis.set_major_formatter(xticks)
ax.legend(ncol=2,fontsize=xy_lb_size)
ax.set_xlabel("Cutoff percentage",fontsize=xy_lb_size)

ax.set_ylabel(ylabels[plt_item],fontsize=xy_lb_size)
plt.savefig(data_path+pic_name,bbox_inches='tight')
plt.show()

# # mag vs radius
# mag = tool_box.mag_generator(5000, 19, 25)
# radii = tool_box.radii_from_mags(mag, 0, 5)
#
# fig = plt.figure(figsize=figs)
# ax = fig.add_subplot(111)
# ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
# for axis in ["bottom", "left", "top", "right"]:
#     ax.spines[axis].set_linewidth(axis_linewidth)
# ax.xaxis.set_tick_params(which="both",direction="in",length=8, width=axis_linewidth)
# ax.yaxis.set_tick_params(which="major",direction="in",length=8, width=axis_linewidth)
# ax.yaxis.set_tick_params(which="minor",direction="in",length=4, width=axis_linewidth)
#
#
# ax.scatter(mag,radii, s=3, color="black")
# ax.set_ylim(0.05, 3)
# x1, x2 = ax.set_xlim(18.9,25.1)
# ax.plot([x1,x2],[0.187, 0.187], color="red",linestyle="--")
#
# ax.set_yscale("log")
# formattor = ScalarFormatter()
# formattor.set_scientific(False)
# ax.set_yticks([0.1, 1])
# ax.set_yticklabels([r"0.1", r"1"])
# ax.yaxis.set_major_formatter(formattor)
# ax.set_xlabel("Magnitude",fontsize=xy_lb_size)
# ax.set_ylabel("Scale length (arcsec)",fontsize=xy_lb_size)
# plt.savefig("E:/works/PICS/mag_radius.png",bbox_inches='tight')
# plt.show()





# # the relation between SNR_F and magnitude
# f = h5py.File("F:/works/selection_bias/data/galsim/dimmerm3/para_0.hdf5")
# mag = f["mag"].value
# f.close()
#
# f = h5py.File("F:/works/selection_bias/data/galsim/dimmerm3/data_0_0.hdf5")
# fsnr = f["data"].value[:,3]
# flux = f["data"].value[:,1]
# mag_m = f["data"].value[:,6]
# idx = flux > 0
# print(idx.sum())
# magm_s = numpy.sort(mag_m[idx])
# numf = idx.sum()
# f.close()
#
# arr = numpy.load("F:/works/selection_bias/data/galsim/dimmerm3/sex_0.npz")["arr_0"]
#
# mag_s = arr[:,2]
# idxs = mag_s > 0
# mags_s = numpy.sort(mag_s[idxs])
# nums = idxs.sum()
#
# fig = plt.figure(figsize=figs)
# ax = fig.add_subplot(111)
# ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
# for axis in ["bottom", "left", "top", "right"]:
#     ax.spines[axis].set_linewidth(axis_linewidth)
# ax.xaxis.set_tick_params(which="both",direction="in",length=4, width=axis_linewidth)
# ax.yaxis.set_tick_params(which="major",direction="in",length=4, width=axis_linewidth)
# ax.yaxis.set_tick_params(which="minor",direction="in",length=2, width=axis_linewidth)
#
#
# num, bins = ax.hist(mag, 100, label="All sample",color="dodgerblue")[:2]
# ys = ax.set_ylim()
# ax.hist(mag[idx],bins, label="Detected by Fourier quad",color="darkorange",alpha=0.8)
# ax.plot([magm_s[int(numf*0.4)],magm_s[int(numf*0.4)]],[ys[0],ys[1]],c="darkorange",linestyle="--")
# ax.hist(mag[idxs],bins, label="Detected by SExtractor",color="limegreen",alpha=0.8)
# ax.plot([mags_s[int(nums*0.4)],mags_s[int(nums*0.4)]],[ys[0],ys[1]],color="limegreen",linestyle="--")
# ax.legend()
# # ax.set_xlim(21,25.1)
# ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
# ax.set_xlabel("Magnitude",fontsize=xy_lb_size)
# ax.set_ylabel("Galaxy number",fontsize=xy_lb_size)
# plt.savefig("F:/works/figs/detection.pdf",bbox_inches='tight')
# plt.show()
