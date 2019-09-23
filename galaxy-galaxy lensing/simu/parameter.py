import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
path.append("D:/Github/astrophy-research/mylib/")
import time
import tool_box
import h5py
import numpy
from plot_tool import Image_Plot
from astropy.cosmology import FlatLambdaCDM


# produce the parameters for GGL simulation
# We use the CMASS w1 sample as the foreground sample now.

# the cosmological parameters
Omega_m0 = 0.31
h = 0.7
H0 = 100*h
cosmos = FlatLambdaCDM(Omega_m0, H0)

# the foreground
h5f = h5py.File("../data/foreground.hdf5","r+")
redshift_f = h5f["/Z"].value
ra_f = h5f["/RA"].value
dec_f = h5f["/DEC"].value
h5f.close()

z_min_f, z_max_f = redshift_f.min(), redshift_f.max()

# the parameter file
h5f = h5py.File("../data/parameter.hdf5", "w")

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

