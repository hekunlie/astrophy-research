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


parent_path = "/mnt/perc/hklee/shear_field/field_1d/"
envs_path = parent_path + "param.dat"

seed = 12356
rng = numpy.random.RandomState(seed)

total_num = 50000
uni_num = 50000

# the background galaxy position
x = numpy.zeros((total_num,))
# symmetrical part
x_uni = rng.uniform(-0, 5, uni_num)
x[:uni_num] = x_uni
# asymmetrical part
if uni_num < total_num:
    x_asy = rng.uniform(1, 4, int(total_num - uni_num))
    x[uni_num:] = x_asy

# shear field
a1, a2 = 0., -0.01
shear_field = a1 + a2*x
print(shear_field.min(), shear_field.max())

# the galactic parameters
h5f = h5py.File(parent_path + "shear_slope.hdf5", "w")

h5f["/param"] = numpy.array([a1,a2])

h5f["/g"] = shear_field

h5f["/x"] = x

mag_i = tool_box.mag_generator(total_num, 20, 24)
flux_i = tool_box.mag_to_flux(mag_i)
h5f["/mag"] = mag_i
h5f["/flux"] = flux_i
h5f.close()
print(flux_i.min(), flux_i.max())


# Plot
img = Image_Plot()
img.subplots(1, 3)
img.set_style()
img.axs[0][0].scatter(x, shear_field)
img.axs[0][1].hist(flux_i, 100)
img.axs[0][2].hist(x, 50)
img.axs[0][0].set_title("shear field")
img.set_label(0, 0, 0, "g")
img.set_label(0, 0, 1, "x")

img.axs[0][1].set_title("Flux hist")
img.axs[0][2].set_title("density")
img.save_img(parent_path + "shear_field_slope.png")
img.show_img()

cmd = "mpirun -np 10 ./simu"
a = Popen(cmd, shell=True)
a.wait()
cmd = "mpirun -np 15 python calculate_chisq.py"
a = Popen(cmd, shell=True)
a.wait()


