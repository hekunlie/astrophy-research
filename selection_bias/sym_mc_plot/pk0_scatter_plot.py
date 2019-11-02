import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/' % my_home)
path.append("D:/Github/astrophy-research/my_lib")
import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
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
ax.yaxis.set_tick_params(which="major",direction="in",length=8, width=axis_linewidth)
ax.yaxis.set_tick_params(which="minor",direction="in",length=4, width=axis_linewidth)




# # the relation between Pk0 and magnitude
#
# f = h5py.File("F:/data_0.hdf5")
# data = f["/data"].value
# detect = data[:, 9]
# select = data[:, 4:7]/64/60
# fsnr = data[:, 0]
# fsnr_f = data[:, 1]
# mag = -data[:, 8]
# f.close()
#
# pk0 = select[:,0]
# pk0_fit = select[:,1]
# # max(pk0, pk0_fit)
# pk0_max = select[:,2]
# idx = detect > -1
#
# mag = mag[idx]
# ch = numpy.random.choice(numpy.array([i for i in range(len(mag))]), 5000, replace=False).tolist()
#
# for axis in ["bottom", "left", "top", "right"]:
#     ax.spines[axis].set_linewidth(axis_linewidth)
# ax.xaxis.set_tick_params(which="both",direction="in",length=7, width=axis_linewidth)
# ax.yaxis.set_tick_params(which="major",direction="in",length=7, width=axis_linewidth)
# ax.yaxis.set_tick_params(which="minor",direction="in",length=4, width=axis_linewidth)
#
# labels = ["P$_{k0}$","P$_{k0}$-fit","MAX(P$_{k0}$,P$_{k0}$-fit)"]
# for i in range(3):
#     ax.scatter(mag[ch], select[:,i][idx][ch], s=5, color="C%d"%i,label=labels[i])
#     # ax.hist(select[:,i],100,label=labels[i])
# ys = ax.set_ylim()
# ax.set_ylim(0.05, ys[1]+300)
# ax.set_yscale("log")
# ax.set_xlabel("Magnitude",fontsize=xy_lb_size)
# ax.set_ylabel("P$_{k0}$",fontsize=xy_lb_size)
# ax.legend(loc="best", bbox_to_anchor=(0.98, 0.98))
# plt.savefig("F:/pk_scatter.pdf",bbox_inches='tight')
# plt.show()


data_path = "G:\galsim\\normal_20_24.8_64x64_galsim_3\sym\sex2_1.5/"
data_path = data_path.replace("\\", "/")
names = ["P$_{k0}$", "P$_{k0}$-fit", "MAX(P$_{k0}$,P$_{k0}$-fit)","MAG$_{true}$", "MAG_AUTO"]
files = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex5", "mag_auto"]

ch_num = 9
cuts_num = 10
x_coord = [i * cuts_num for i in range(ch_num)]
ch = [i for i in range(ch_num)]

plt_item = 0

ylabels = ["m$_1 \\times 10^2$", "m$_2 \\times 10^2$",
           "c$_1 \\times 10^4$", "c$_2 \\times 10^4$"]
mc_item = ["m1", "m2", "c1", "c2"]
pic_name = data_path + "pk_comp_%s.pdf"%mc_item[plt_item]
for i in range(len(files)):
    npz_path = data_path+files[i]+"/total.npz"
    print(os.path.exists(npz_path))
    data = numpy.load(npz_path)
    mc1 = data['arr_0'][:, ch]
    mc2 = data['arr_1'][:, ch]
    if i == 2:
        line_style = "-"
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
ax.set_ylim(ys[0], ys[1]+0.5)
ax.xaxis.set_major_formatter(xticks)
ax.legend(ncol=2,fontsize=xy_lb_size)
ax.set_xlabel("Cutoff percentage",fontsize=xy_lb_size)

ax.set_ylabel(ylabels[plt_item],fontsize=xy_lb_size)
plt.savefig(pic_name,bbox_inches='tight')
plt.show()