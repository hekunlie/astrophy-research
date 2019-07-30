import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
import tool_box
import h5py
import time

# set up Z bins for foreground source
# foreground name
# area name: w_1 ..
# Z bin: [z1, z2]
# new file name: w_1_1 ..

fore_source = argv[1]
area = argv[2]
z1, z2 = float(argv[3]), float(argv[4])
new_name = argv[5]
parent_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/foreground/%s/"%fore_source
data_path = "%s/%s.hdf5"%(parent_path, area)
pic_path = "%s/%s.png"%(parent_path, new_name)

h5f = h5py.File(data_path,"r")
names = list(h5f.keys())
redshift = h5f["/Z"].value
num = redshift.shape[0]

idx1 = redshift >= z1
idx2 = redshift < z2
idx = idx1 & idx2

img = Image_Plot()
img.subplots(1,1)
bins = img.axs[0][0].hist(redshift, label="All")[1]
img.axs[0][0].hist(redshift[idx], bins=bins, label="[%.3f, %.3f]"%( z1, z2))
img.save_img(pic_path)

h5f_new = h5py.File(parent_path + new_name + ".hdf5","w")
for nm in names:
    h5f_new["/%s"%nm] = h5f["/%s"%nm].value[idx]
    print("%d == > %d in [%.3f, %.3f] Z bin"%(num, idx.sum(), z1, z2))
h5f.close()
h5f_new.close()

