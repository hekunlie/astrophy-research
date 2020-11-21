import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
path.append("D:/GitHub/astrophy-research/mylib")
# from plot_tool import Image_Plot
import h5py
import numpy


tag = argv[1]
h5f = h5py.File("/home/hklee/work/CFHT/correlation/result/core_%s_num_count.hdf5"%tag,"r")
data = h5f["/0/data"][()]
jack = h5f["/0/jack_label"][()]
h5f.close()
block_num = int(jack.shape[0]/2)
print(block_num)

pts_num = 147

theta = data[:,:pts_num].sum(axis=0)
num = data[:,pts_num:int(pts_num*2)].sum(axis=0)
theta_mean = (theta/num).reshape(21,7)
print(theta_mean)
print(theta_mean.shape)

# img = Image_Plot()
# img.subplots(1,1)
# for i in range(21):
#     img.axs[0][0].scatter(theta_mean[i],range(7))
# img.axs[0][0].set_xscale("log")
# img.show_img()