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



# all the coordinates should be converted to the unit of arcmin
# number of grid
nx = 200
ny = nx
# the field size in unit of arcmin
delta_ra = 20
delta_dec = delta_ra
half_side = nx/2
# arcmin/pix
pixel_scale = delta_ra/nx
# the shear profile
sigma = 2 # arcmin
amplitude = 0.8
dx, dy = 0, 0
# galaxy number density
dens_num = 75
num_each_expo = dens_num*delta_ra*delta_dec
expo_num = 50
total_num = num_each_expo*expo_num
profile_params = [amplitude, dx, dy, sigma]
print("%d galaxies each exposure. %d exposures."%(num_each_expo, expo_num))

parent_path = "/mnt/perc/hklee/CFHT/multi_shear/cluster_field/"

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

# shear field
g = MCMC_program.shear_field(profile_params, x, y)


img = Image_Plot()
img.subplots(2,2)
inverse = range(nx-1, -1, -1)
for i in range(2):
    for j in range(2):
        tag = i*2 + j
        if tag < 3:
            fig = img.axs[i][j].imshow(g[tag][inverse],cmap="jet")
            img.figure.colorbar(fig, ax=img.axs[i][j])
        else:
            num_dens = numpy.histogram2d(gal_coord[0], gal_coord[1], [ra_bin, dec_bin])[0]
            fig = img.axs[i][j].imshow(num_dens,cmap="jet")
            img.figure.colorbar(fig, ax=img.axs[i][j])

        img.set_label(i, j, 0, "DEC. [arcmin]")
        img.set_label(i, j, 1, "R.A. [arcmin]")
        img.axs[i][j].set_yticks([int(nx - k * nx / 5) for k in range(6)])
        img.axs[i][j].set_yticklabels(["%.1f" % (k * nx / 5. * pixel_scale) for k in range(6)])
        img.axs[i][j].set_xticks([int(k * ny / 5) for k in range(6)])
        img.axs[i][j].set_xticklabels(["%.1f" % (k * ny / 5. * pixel_scale) for k in range(6)])
img.save_img(parent_path + "pic/shear_field.png")
img.show_img()

# the galactic parameters
h5f = h5py.File(parent_path + "param.hdf5", "w")

h5f["/FIELD_g"] = g[0]
h5f["/FIELD_g1"] = g[1]
h5f["/FIELD_g2"] = g[2]

ra = gal_coord[0]
dec = gal_coord[1]
print(ra.min(), ra.max(), ra.mean())
print(dec.min(), dec.max(), dec.mean())

h5f["/ra"] = ra
h5f["/dec"] = dec
h5f["/radius"] = gal_coord[2]

g, g1, g2 = MCMC_program.shear_field(profile_params, ra, dec)
h5f["/g"] = g
h5f["/g1"] = g1
h5f["/g2"] = g2
print(g1.min(), g1.max())

mag_i = tool_box.mag_generator(num_each_expo, 20, 23.5)
flux_i = tool_box.mag_to_flux(mag_i)
h5f["/mag"] = mag_i
h5f["/flux"] = flux_i
print(flux_i.min(), flux_i.max())

img = Image_Plot()
img.subplots(1, 3)

norm_g = plt.Normalize(vmin=numpy.min(g), vmax=numpy.max(g))
cmap_g = plt.get_cmap('jet')
img.axs[0][0].scatter(ra, dec, color=cmap_g(norm_g(g)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][0])

norm_g = plt.Normalize(vmin=numpy.min(g1), vmax=numpy.max(g1))
cmap_g = plt.get_cmap('jet')
img.axs[0][1].scatter(ra, dec, color=cmap_g(norm_g(g1)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][1])

norm_g = plt.Normalize(vmin=numpy.min(g2), vmax=numpy.max(g2))
cmap_g = plt.get_cmap('jet')
img.axs[0][2].scatter(ra, dec, color=cmap_g(norm_g(g2)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][2])

img.save_img(parent_path + "pic/expo.png")
img.close_img()
h5f.close()





