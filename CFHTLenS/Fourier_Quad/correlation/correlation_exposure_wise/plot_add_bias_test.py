import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
# path.append('%s/work/mylib/' % my_home)
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import h5py
import numpy
import matplotlib.pyplot as plt
import tool_box


data_path = "E:/works/correlation/CFHT/2021_6_15/cut_2.5_new/unsmooth/ODDS_0.4"

h5f = h5py.File(data_path + "/bias_test.hdf5","r")

img = Image_Plot(fig_x=4,fig_y=3,xpad=0.14,ypad=0.2)
img.subplots(2,4)

for i in range(2):
    for j in range(2):
        tag = int(2*i + j)
        data = numpy.abs(h5f["/%d"%tag][()])
        scale = 1
        # img.axs[0][tag].errorbar(data[:,0], data[:,1]*scale, data[:,2]*scale,fmt=" ",capsize=3,marker="s", label="$10^3 c_1$")
        img.axs[0][tag].scatter(data[:,0], data[:,1]*scale/data[:,2]*scale,marker="s", label="$|c_1|/\delta_{c_1}$")
        # img.axs[1][tag].errorbar(data[:,0], data[:,3]*scale, data[:,4]*scale,fmt=" ",capsize=3,marker="s", label="$10^3 c_2$")
        img.axs[1][tag].scatter(data[:,0], data[:,3]*scale/data[:,4]*scale,marker="s", label="$|c_2|/\delta_{c_2}$")
        img.set_label(0,tag,1,"Jack_id")
        img.set_label(1,tag,1,"Jack_id")
        img.axs[0][tag].legend()
        img.axs[1][tag].legend()
        # img.axs[0][tag].set_yscale("log")
        # img.axs[1][tag].set_yscale("log")
h5f.close()

img.save_img(data_path + "/add_bias_test.jpg")
img.show_img()