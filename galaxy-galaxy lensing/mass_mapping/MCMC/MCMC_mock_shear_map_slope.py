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



parent_path = "/mnt/perc/hklee/CFHT/multi_shear/cluster_field/"
envs_path = "%s/param_slope.dat" % my_home

# all the coordinates should be converted to the unit of arcmin
# number of grid
nx = 400
ny = nx
# the field size in unit of arcmin
delta_ra = 16
delta_dec = delta_ra
half_side = nx/2
# arcmin/pix
pixel_scale = delta_ra/nx
# galaxy number density
dens_num = 75
num_each_expo = dens_num*delta_ra*delta_dec
expo_num = 30
total_num = num_each_expo*expo_num
print("%d galaxies each exposure. %d exposures."%(num_each_expo, expo_num))

# the grid
# x: ra, y: dec
my, mx = numpy.mgrid[-half_side:half_side,-half_side:half_side]
y, x = my*pixel_scale, mx*pixel_scale

ra_bin = numpy.linspace(-delta_ra/2, delta_ra/2, nx+1)
dec_bin = numpy.linspace(-delta_dec/2, delta_dec/2, ny+1)

ra_range = (-delta_ra/2, delta_ra/2)
dec_range = (-delta_dec/2, delta_dec/2)

# the background galaxy position
gal_coord = numpy.zeros((3, num_each_expo))
gal_coord[0] = numpy.random.uniform(ra_range[0], ra_range[1], num_each_expo)
gal_coord[1] = numpy.random.uniform(dec_range[0], dec_range[1], num_each_expo)
gal_coord[2] = numpy.sqrt(gal_coord[0]**2 + gal_coord[1]**2)


a1, a2, a3 = 0, 0.002, -0.006
profile_params = [a1, a2, a3]
# shear field
shear_field = MCMC_program.shear_slope(profile_params, x, y)

contents = [['param', "grid_nx", '%d'%nx], ['param', "grid_ny", '%d'%ny],
            ['param', "RA", '%.2f'%delta_ra], ['param',"DEC", '%.2f'%delta_dec],
            ['param', "pixel_scale", '%.6f' % pixel_scale], ['param', "a1", '%.4f' % a1],
            ['param', "a2", '%.2f' % a2], ['param', "a3", '%.4f' % a3],
            ['param', "density/arcmin^2", '%d' %dens_num], ['param', "exposure", '%d'%expo_num]]

path_items = tool_box.config(envs_path, ['add' for i in range(len(contents))], contents, True)


# the galactic parameters
h5f = h5py.File(parent_path + "param_slope.hdf5", "w")

h5f["/FIELD_g"] = shear_field

ra = gal_coord[0]
dec = gal_coord[1]
print(ra.min(), ra.max(), ra.mean())
print(dec.min(), dec.max(), dec.mean())

h5f["/ra"] = ra
h5f["/dec"] = dec
h5f["/radius"] = gal_coord[2]

shear_field_ch = MCMC_program.shear_slope(profile_params, ra, dec)
h5f["/g"] = shear_field_ch
print(shear_field_ch.min(), shear_field_ch.max())

mag_i = tool_box.mag_generator(num_each_expo, 20, 23.5)
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





