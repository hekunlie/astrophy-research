import matplotlib
import matplotlib.pyplot as plt
# matplotlib.use("Agg")
import os
from sys import path
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
# path.append('%s/work/mylib/'%my_home)
path.append('D:/GitHub/astrophy-research/mylib')
import numpy
from plot_tool import Image_Plot
import h5py

h5f = h5py.File("E:/mag_true_12.hdf5", "r")
mag_true = h5f["/data"][()]
h5f.close()

h5f = h5py.File("E:/mask_12.hdf5", "r")
mask = h5f["/data"][()]
h5f.close()
idx = mask > 0
print(idx.sum()/idx.shape[0])

matplotlib.style.use('default')
plt.rcParams['font.family'] = 'serif'

scale = [24.2423, 23.4742]
line_labels = ["$30%$ cutoff","$60%$ cutoff"]
img = Image_Plot(legend_size=14, fig_x=6, fig_y=4,plt_line_width=2,axis_linewidth=2,xpad=0.2)
img.subplots(1,1)
# img.set_style()
img.axis_type(0,"major",tick_len=6, tick_width=1.2)
img.axis_type(1,"major",tick_len=6, tick_width=1.2)
bins = img.axs[0][0].hist(-mag_true, 50, ec="gray", label="Entire sample")[1]
img.axs[0][0].hist(-mag_true[idx], bins, ec="gray", label="Detected sample")

ys = img.axs[0][0].set_ylim()
for i,s in enumerate(scale):
    img.axs[0][0].plot([s,s],[ys[0], ys[1]], lw=img.plt_line_width, ls="--", label=line_labels[i])

img.axs[0][0].legend(fontsize=img.legend_size,frameon=False)
img.set_ticklabel_str(0,0,0,[2*i*100000 for i in range(1,7)],["%d"%(2*i) for i in range(1,7)])
img.set_label(0,0,0,"$10^{-5}N$")
img.set_label(0,0,1,"Magnitude")
# img.axs[0][0].set_xlim(22, 25)
img.save_img("E:/cutoff_line.pdf")
img.show_img()
