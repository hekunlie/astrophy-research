import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import numpy
import h5py
from plot_tool import Image_Plot
import tool_box
import matplotlib.pyplot as plt


size = 30
sig = 5
cen = int(size/2)
cut = 6
shear = tool_box.gauss_profile(size, sig, cen, cen)*9 - 0.025
g1 = shear[cen-cut:cen,cen-cut:cen]
g2 = shear[cen-cut:cen,cen:cen+cut]
g = numpy.sqrt(g1**2+g2**2)
print(g1.std(),g2.std())

h5f = h5py.File("/mnt/perc/hklee/CFHT/multi_shear/data/shear.hdf5","w")
h5f["/g1"] = g1
h5f["/g2"] = g2
h5f.close()

img = Image_Plot()
img.subplots(1, 2)
img.set_style()
cmap = plt.get_cmap('jet')
fig = img.axs[0][0].imshow(g1, cmap=cmap)
img.axs[0][0].set_title("STD: %.4f"%g1.std())
img.figure.colorbar(fig, ax=img.axs[0][0])

img.axs[0][1].set_title("STD: %.4f"%g2.std())
fig = img.axs[0][1].imshow(g2, cmap=cmap)
img.figure.colorbar(fig, ax=img.axs[0][1])
img.save_img("/mnt/perc/hklee/CFHT/multi_shear/data/shear.png")
# img.show_img()