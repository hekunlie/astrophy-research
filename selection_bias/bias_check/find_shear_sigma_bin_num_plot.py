import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box

img = Image_Plot(fig_x=6, fig_y=8,xpad=0.1,ypad=0.2)
img.subplots(1,3)
signal = 0.01
bin_num = [2,4,8, 12, 14, 16,20,32,64, 96,128,192]
data_num = len(bin_num)
y = range(data_num)
titles = ["Num = $10^5$,$\sigma=1$,$\mu=0.01$","Num = $10^6$,$\sigma=1$,$\mu=0.01$","Num = $10^7$,$\sigma=1$,$\mu=0.01$"]
scale = [100,100,1000]
x_ticks = numpy.linspace(-0.02,0.01,4)
for i in range(3):
    h5f = h5py.File("D:/data_%d.hdf5"%(i+1),"r")
    data = h5f["/data"][()]
    img.axs[0][i].errorbar(data[1,:data_num],y,xerr=data[2,:data_num],capsize=img.cap_size-1)
    ys = img.axs[0][i].set_ylim()
    img.axs[0][i].plot([signal,signal],[ys[0],ys[1]], ls="dotted", c="grey", alpha=0.5)
    img.set_label(0, i, 1, "g")
    for j in range(data_num):
        if i < 2:
            text_str = "%d bins, $10^2\sigma=%.3f$, $N\sigma^2=%.3f$"%(data[4,j],scale[i]*data[2,j],data[3,j])
        else:
            text_str = "%d bins, $10^3\sigma=%.3f$, $N\sigma^2=%.3f$"%(data[4,j],scale[i]*data[2,j],data[3,j])
        img.axs_text(0, i, y[j], -0.019, text_str,text_fontsize=img.legend_size-4,ax_trans=False)
    img.axs[0][i].set_xlim(signal-0.031,signal+0.007)
    img.axs[0][i].set_title(titles[i],fontsize=img.legend_size)
    img.del_tick(0,i,[0])
    img.axs[0][i].set_xticks(x_ticks)
img.save_img("D:/sigma.png")
img.show_img()