import matplotlib
# matplotlib.use("Agg")
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import numpy
import h5py


parent_path = "E:/works/Group_meeting/2019-11-25-shear_bias_checking"
img = Image_Plot(fig_x=6,fig_y=4,ypad=0.2)
img.subplots(1,2)
img.set_style()
for i in range(6):
    h5f = h5py.File(parent_path+"/result_%d.hdf5"%i,"r")
    mc = h5f["/mc"].value
    data = h5f["/data"].value
    h5f.close()
    tag = 2.5*i+2
    img.axs[0][0].errorbar(mc[0,2]-1, tag, xerr=mc[0,3],c="C%d"%i,
                           marker="v",mfc="none",capsize=img.cap_size)
    img.axs[0][0].errorbar(mc[0,6]-1, tag+0.4, xerr=mc[0,7],c="C%d"%i,
                           marker="o",mfc="none",capsize=img.cap_size)

    img.axs[0][0].errorbar(mc[1,2]-1, tag, xerr=mc[1,3],c="C%d"%i,
                           marker="v",capsize=img.cap_size)
    img.axs[0][0].errorbar(mc[1,6]-1, tag+0.4, xerr=mc[1,7],c="C%d"%i,
                           marker="o",capsize=img.cap_size)

    img.axs[0][1].errorbar(10000*mc[0,0], tag, xerr=10000*mc[0,1],c="C%d"%i,
                           marker="v",mfc="none",capsize=img.cap_size)
    img.axs[0][1].errorbar(10000*mc[0,4], tag+0.45, xerr=10000*mc[0,5],c="C%d"%i,
                           marker="o",mfc="none",capsize=img.cap_size)

    img.axs[0][1].errorbar(10000*mc[1,0], tag+0.1, xerr=10000*mc[1,1],c="C%d"%i,
                           marker="v",capsize=img.cap_size)
    img.axs[0][1].errorbar(10000*mc[1,4], tag+0.55, xerr=10000*mc[1,5],c="C%d"%i,
                           marker="o",capsize=img.cap_size)
for i in range(2):
    img.axs[0][i].plot([0, 0], [-2, 16], c="grey", ls="dashed")
    img.axs[0][i].set_ylim(0, 14)

img.axs[0][0].set_xlim(-0.0063,0.0063)
img.axs[0][1].set_xlim(-0.92,0.92)


img.set_label(0,0,1,"$m$")
img.set_label(0,1,1,"$10^4 c$")
img.save_img(parent_path+"/m_scatter.png")
img.show_img()