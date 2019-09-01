import matplotlib
matplotlib.use("Agg")
import numpy
from sys import path
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
import h5py
from subprocess import Popen
from Fourier_Quad import Fourier_Quad


parent_path = "/mnt/perc/hklee/multi_shear_dectect/"
shear_num = 10
source_num = 800000
for i in range(shear_num):

    h5f = h5py.File(parent_path + "param_%d.hdf5"%i, "w")
    mag = tool_box.mag_generator(source_num, 20, 24)
    flux = tool_box.mag_to_flux(mag)
    h5f["/flux"] = flux
    h5f.close()

print("Begin simu...")
cmd = "mpirun -np 20 ./simu"
a = Popen(cmd, shell=True)
a.wait()

fq = Fourier_Quad(3, 123)

for i in range(shear_num):
    h5f = h5py.File(parent_path + "result_%d.hdf5"%i, "r")
    data = h5f["/data"].value
    mg1 = data[:,0]
    mg2 = data[:,1]
    mnu1 = data[:,2] + data[:,3]
    mnu2 = data[:,2] - data[:,3]
    h5f.close()

    gh1, gh1_sig = fq.find_shear_new(mg1,mnu1, 8)[:2]
    gh2, gh2_sig = fq.find_shear_new(mg2,mnu2, 8)[:2]

    g1 = -0.05 + i * 0.01
    g2 = 0.05 - i * 0.01
    print("%.5f: %.5f (%.5f), %.5f: %.5f (%.5f)"%(g1, gh1, gh1_sig, g2, gh2, gh2_sig))