import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/'%my_home)
import numpy
import tool_box
from mpi4py import MPI
import h5py
from plot_tool import Image_Plot
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


data_path = argv[1]

chip_num = 1500
stamp_num = 10000

t1 = time.time()

# magnitude & flux
mag_s, mag_e = 20.5, 25.5
mag = tool_box.mag_generator(chip_num*stamp_num, mag_s, mag_e).astype(dtype=numpy.float32)
flux = tool_box.mag_to_flux(mag).astype(dtype=numpy.float32)

# galactic radius
radius_s, radius_e = 0.75, 1.87

radius = tool_box.radii_from_mags(mag, radius_s, radius_e).astype(dtype=numpy.float32)/0.187

if rank == 0:
    if not os.path.exists(data_path):
        os.makedirs(data_path)
    if not os.path.exists(data_path + "/pic"):
        os.makedirs(data_path + "/pic")

    g = numpy.random.uniform(0, 0.04, 15)
    theta = numpy.random.uniform(0, numpy.pi, 15)
    g1 = g * numpy.cos(2 * theta)
    g2 = g * numpy.sin(2 * theta)
    h5f = h5py.File(data_path + "/shear.hdf5", "w")
    h5f["/g1"] = g1
    h5f["/g2"] = g2
    h5f.close()

    img = Image_Plot()
    img.subplots(1, 1)
    img.axs[0][0].scatter(g1, g2,c="k")
    img.save_img(data_path + "/pic/shear.png")
    img.close_img()

comm.Barrier()

h5f = h5py.File(data_path + "/paras_%d.hdf5"%rank,"w")
h5f["/flux"] = flux
h5f["/mag"] = mag
h5f["/radius"] = radius
h5f.close()

img = Image_Plot()
img.subplots(1, 2)
img.axs[0][0].hist(mag, 100, histtype="step")
img.axs[0][1].scatter(mag[:5000], radius[:5000],c="k",s=5)
img.save_img(data_path + "/pic/para_%d.png"%rank)
img.close_img()

comm.Barrier()
t2 = time.time()
if rank == 0:
    print(t2-t1)

