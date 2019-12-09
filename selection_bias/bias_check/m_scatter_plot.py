import matplotlib
# matplotlib.use("Agg")
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import numpy
import h5py


parent_path = "E:/works/Group_meeting/2019-11-25-shear_bias_checking/result/data_from_pi/3"
img = Image_Plot(fig_x=6,fig_y=8, ypad=0.1)
img.subplots(1,2)
img.set_style()

color_tag = [0]

x = numpy.linspace(-1, 1, 100)
pts_size = 7
for i in range(12):
    if i>0:
        h5f = h5py.File(parent_path+"/%d/result_noisy_bin_num_40.hdf5"%(i-1),"r")
    else:
        h5f = h5py.File(parent_path + "/result_all/result_noisy.hdf5", "r")
    mc = h5f["/sym_mc"].value
    h5f.close()

    tag = 3*i + 2
    if i in color_tag:
        color = "C3"
    else:
        color = "C2"

    # m
    img.axs[0][0].fill_between(x, tag - 0.8, tag + 1.2, facecolor='grey', alpha=0.3)

    img.axs[0][0].errorbar(mc[0,0], tag, xerr=mc[0,1],c=color,
                           marker="v",mfc="none",capsize=img.cap_size,ms=pts_size)
    img.axs[0][0].errorbar(mc[1,0]-1, tag+0.6, xerr=mc[1,1],c=color,
                           marker="o",mfc="none",capsize=img.cap_size,ms=pts_size)

    img.axs[0][0].errorbar(mc[1,2]-1, tag, xerr=mc[1,3],c=color,
                           marker="v",capsize=img.cap_size,ms=pts_size)
    img.axs[0][0].errorbar(mc[1,6]-1, tag+0.6, xerr=mc[1,7],c=color,
                           marker="o",capsize=img.cap_size,ms=pts_size)


    # c
    img.axs[0][1].fill_between(x, tag-0.8, tag + 1.2,facecolor='grey', alpha=0.3)
    img.axs[0][1].errorbar(10000*mc[0,0], tag, xerr=10000*mc[0,1],c=color,
                           marker="v",mfc="none",capsize=img.cap_size,ms=pts_size)
    img.axs[0][1].errorbar(10000*mc[0,4], tag+0.6, xerr=10000*mc[0,5],c=color,
                           marker="o",mfc="none",capsize=img.cap_size,ms=pts_size)

    img.axs[0][1].errorbar(10000*mc[1,0], tag, xerr=10000*mc[1,1],c=color,
                           marker="v",capsize=img.cap_size,ms=pts_size)
    img.axs[0][1].errorbar(10000*mc[1,4], tag+0.6, xerr=10000*mc[1,5],c=color,
                           marker="o",capsize=img.cap_size,ms=pts_size)


for i in range(2):
    img.axs[0][i].plot([0, 0], [-2, tag+6], c="grey", ls="dashed")
    if i == 0:
        img.axs[0][i].plot([0.002, 0.002], [-2, tag+6], c="grey", ls="dashed",alpha=0.5)
        img.axs[0][i].plot([-0.002, -0.002], [-2, tag+6], c="grey", ls="dashed",alpha=0.5)
    img.axs[0][i].set_ylim(0, tag+3)
    img.axs[0][i].set_yticklabels([])

img.axs[0][0].set_xlim(-0.008,0.0061)
img.axs[0][1].set_xlim(-0.92,0.92)


img.set_label(0,0,1,"$m$")
img.set_label(0,1,1,"$10^4 c$")
img.save_img(parent_path+"/m_scatter.png")
img.show_img()