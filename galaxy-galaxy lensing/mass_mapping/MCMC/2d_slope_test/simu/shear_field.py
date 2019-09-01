import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
from sys import path,argv
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
import h5py
import MCMC_program


parent_path = "/mnt/perc/hklee/shear_field/field_1d/"
envs_path = parent_path + "param.dat"

total_num = 20000
x_st, x_ed = -8, 8
# the background galaxy position
x = numpy.random.uniform(x_st, x_ed, total_num)

a1, a2 = 0, -0.02
profile_params = [a1, a2]

# shear field
shear_field = a1 + a2*x

# the galactic parameters
h5f = h5py.File(parent_path + "param_slope.hdf5", "w")

h5f["/FIELD_g"] = shear_field

h5f["/x"] = x


mag_i = tool_box.mag_generator(total_num, 20, 23.5)
flux_i = tool_box.mag_to_flux(mag_i)
h5f["/mag"] = mag_i
h5f["/flux"] = flux_i
print(flux_i.min(), flux_i.max())

# Plot
img = Image_Plot()
img.subplots(1, 3)
img.set_style()
inverse = range(nx - 1, -1, -1)

norm_g = plt.Normalize(vmin=numpy.min(shear_field_ch), vmax=numpy.max(shear_field_ch))
cmap_g = plt.get_cmap('jet')

for i in range(3):
    if i == 0:
        fig = img.axs[0][i].imshow(shear_field[inverse], cmap="jet")
        img.figure.colorbar(fig, ax=img.axs[0][i])
    elif i == 1:
        num_dens = numpy.histogram2d(gal_coord[0], gal_coord[1], [ra_bin, dec_bin])[0]
        fig = img.axs[0][i].imshow(num_dens, cmap="jet")
        img.figure.colorbar(fig, ax=img.axs[0][i])
    else:
        img.axs[0][i].scatter(ra, dec, color=cmap_g(norm_g(shear_field_ch)),s=1)
        sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
        sm._A = []
        plt.colorbar(sm, ax=img.axs[0][i])
        img.axs[0][i].set_title("%d galaxies"%num_each_expo)

    img.set_label(0, i, 0, "DEC. [arcmin]")
    img.set_label(0, i, 1, "R.A. [arcmin]")
    if i < 2:
        img.axs[0][i].set_yticks([int(nx - k * nx / 5) for k in range(6)])
        img.axs[0][i].set_yticklabels(["%.1f" % (k * nx / 5. * pixel_scale) for k in range(6)])
        img.axs[0][i].set_xticks([int(k * ny / 5) for k in range(6)])
        img.axs[0][i].set_xticklabels(["%.1f" % (k * ny / 5. * pixel_scale) for k in range(6)])

img.save_img(parent_path + "pic/shear_field_slope.png")
img.show_img()





