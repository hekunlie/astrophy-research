import numpy
import matplotlib.pyplot as plt
from sys import path
# import os
# my_home = os.popen("echo $HOME").readlines()[0][:-1]
# path.append('%s/work/mylib/'%my_home)
path.append("D:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box


parent_path = "G:/cata_result/cfht/fore_bin_1/"

img = Image_Plot()
img.subplots(1,1)

file_num = 6
bin_num = 13

x = tool_box.set_bin_log(0.04, 15, bin_num+1)
print(x)
signal = numpy.zeros((file_num,bin_num))
sigma = numpy.zeros((file_num,bin_num))
los_dist = numpy.zeros((file_num,bin_num))

for i in range(file_num):
    data_path = parent_path + "cmass_result_w_1_%s.npz"%i
    data = numpy.load(data_path)["arr_0"]
    signal[i] = data[0]
    sigma[i] = data[1]
    los_dist[i] = data[-1]

for i in range(6):
    print(x.shape, signal[:,i].shape)
    img.axs[0][0].errorbar(los_dist[:,i], signal[:,i], sigma[:,i], label="[%.5f, %.5f] Mpc/h"%(x[i], x[i+1]), capsize=3)
img.axs[0][0].legend()
img.set_label(0,0,0,"$\Delta\Sigma\  [\\rm{h\cdot M_{\odot}\cdot pc^{-2}}]$",size=15)
img.set_label(0,0,1,"LOS distance$\ [\\rm{Mpc \cdot h^{-1}}]$",size=15)
img.show_img()
