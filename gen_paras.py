import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
from Fourier_Quad import Fourier_Quad
import tool_box
import h5py


with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "=" in path:
        env_location, env_path = path.split("=")[0:2]
        if "cf_data_path" == env_location:
            cf_data_path = env_path

with open("./paras.dat","r") as f:
    contents = f.readlines()
for item in contents:
    if ":" in item:
        item_name, item_val = item.split(":")
        if "stamp_size" == item_name:
            stamp_size = int(item_val)
        if "total_num" == item_name:
            num = int(item_val)
        if "g1_sig" == item_name:
            g1_sig = float(item_val)
        if "g2_sig" == item_name:
            g2_sig = float(item_val)
        if "g1_correlation" == item_name:
            gg_cor1 = float(item_val)
        if "g2_correlation" == item_name:
            gg_cor2 = float(item_val)
        if "g1_range_s" == item_name:
            g1_s = float(item_val)
        if "g1_range_e" == item_name:
            g1_e = float(item_val)
        if "g2_range_s" == item_name:
            g2_s = float(item_val)
        if "g2_range_e" == item_name:
            g2_e = float(item_val)
        if "magnitude_s" == item_name:
            mag_s = int(item_val)
        if "magnitude_e" == item_name:
            mag_e = int(item_val)

para_path = cf_data_path + "para.hdf5"
f = h5py.File(para_path,'w')

# generate the correlated (g1, g2) pairs

g1_pairs = tool_box.rand_gauss2([g1_s, g1_e], [g1_s, g1_e], num, g1_sig, g1_sig, gg_cor1)
g2_pairs = tool_box.rand_gauss2([g2_s, g2_e], [g2_s, g2_e], num, g2_sig, g2_sig, gg_cor2)

plt.figure(figsize=(10,5))
plt.subplot(121)
plt.hist2d(g1_pairs[:,0],g1_pairs[:,1], 50)
plt.colorbar()
plt.subplot(122)
plt.hist2d(g2_pairs[:,0],g2_pairs[:,1], 50)
plt.colorbar()
gg_cor_fig = cf_data_path + "gg_cor.png"
plt.savefig(gg_cor_fig)
plt.close()

g = numpy.row_stack((g1_pairs, g2_pairs))
f["/shear"] = g

# ellipticity

ellips = numpy.zeros((2*num, 2))
for i in range(2):
    seed = 123 + i
    rng = numpy.random.RandomState(seed=seed)
    e = tool_box.ellip_mock(num, seed)
    theta = rng.uniform(0, numpy.pi, num)
    q = (1 - e) / (1 + e)
    es = (1 - q ** 2) / (1 + q ** 2)
    e1 = es * numpy.cos(2 * theta)
    e2 = es * numpy.sin(2 * theta)
    ellips[i*num: (i+1)*num, 0] = e1
    ellips[i*num: (i+1)*num, 1] = e2

    plt.subplot(131)
    plt.hist(e, 50)
    plt.subplot(132)
    plt.hist(e1, 50)
    plt.subplot(133)
    plt.hist(e2, 50)
    pic_name = cf_data_path + "ellip_%d.png"%i
    plt.savefig(pic_name)
    plt.close()

f["/ellips"] = ellips

# magnitude
magnitude = numpy.zeros((2*num, 1))
for i in range(2):
    mags = tool_box.mags_mock(num, mag_s, mag_e)
    magnitude[i*num: (i+1)*num, 0] = mags

f["/magnitudes"] = magnitude

f.close()