import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
path.append("D:/Github/astrophy-research/mylib/")
import tool_box
import h5py
import numpy
from plot_tool import Image_Plot

from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/single_shear_test"

# the parameter file
h5f = h5py.File(data_path + "/sheared_para_.hdf5", "w")

# radius bin
radius_bin_num = 15
radius_bin = tool_box.set_bin_log(0.04, 15, radius_bin_num)
h5f["/radius_bin"] = radius_bin

# source num in each radius bin
num_each_bin = numpy.zeros((radius_bin_num-1,), dtype=numpy.intc)
for i in range(radius_bin_num - 1):
    num_each_bin[i] = int(20000*i*i + 10000)
h5f["num_each_bin"] = num_each_bin

#

print(num_each_bin)
img = Image_Plot()
img.subplots(1,1)
img.axs[0][0].plot(radius_bin[1:], num_each_bin)
img.show_img()

