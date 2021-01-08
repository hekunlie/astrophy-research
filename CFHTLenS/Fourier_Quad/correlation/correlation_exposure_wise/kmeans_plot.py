import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from plot_tool import Image_Plot


src_path = argv[1]

ra_bin = [1000,5000,10000,15000,30000]
labels = []
sub_num = numpy.array([75,34,56,35],dtype=numpy.intc)
img = Image_Plot(fig_x=8, fig_y=6)
img.subplots(2,2)
for i in range(2):
    for j in range(2):
        tag = int(i*2 + j)

        h5f = h5py.File(src_path + "/group_predict_%d.hdf5"%tag, "r")
        group_label = h5f["/data"][()]
        ra_dec = h5f["/ra_dec"][()]
        h5f.close()

        for k in range(sub_num[tag]):
            idx = group_label == k
            print(k, idx.sum())
        # group_label += sub_num[:tag].sum()
        print(tag, group_label.min(), group_label.max(),sub_num[:tag].sum())
        print("\n")
        img.pts_scatter(i,j,ra_dec[:,0],ra_dec[:,1],group_label,[0,sub_num[tag]])
img.save_img("kmeans_new.png")
img.close_img()

