import matplotlib
matplotlib.use("Agg")
import numpy
from sys import path,argv
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
import h5py


# all the coordinates should be converted to the unit of arcmin
def shear_field(x, y, cen, sig, ampl):
    ''' the shear field '''
    r = numpy.sqrt(x ** 2 + y ** 2)
    g = ampl*r*numpy.exp(-((y - shift[0]) ** 2 + (x-shift[1]) ** 2) / 2 / sig / sig) / numpy.pi / 2 / sig/sig
    if type(x) is not float:
        idx = numpy.abs(g) < 0.001
        g[idx] = 0
        idx = r == 0
        r[idx] = 1
        x[idx] = 0
        y[idx] = 0
    sin_theta = x/r
    cos_theta = y/r
    sin_2theta = 2*sin_theta*cos_theta
    cos_2theta = cos_theta**2 - sin_theta**2
    g1 = g*cos_2theta
    g2 = -g*sin_2theta
    return g, g1, g2


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
sigma = 1.5 # arcmin
amplitude = 0.8
shift = (0, 0)
# galaxy number density
dens_num = 25
num_each_expo = dens_num*delta_ra*delta_dec
expo_num = 6
total_num = num_each_expo*expo_num
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
gal_coord = numpy.zeros((3, total_num))
gal_coord[0] = numpy.random.uniform(ra_range[0], ra_range[1], total_num)
gal_coord[1] = numpy.random.uniform(dec_range[0], dec_range[1], total_num)
gal_coord[2] = numpy.sqrt(gal_coord[0]**2 + gal_coord[1]**2)

# shear field
g = shear_field(x, y, shift, sigma, amplitude)


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
            print(num_dens.sum(), total_num)

        img.set_label(i, j, 0, "DEC. [arcmin]")
        img.set_label(i, j, 1, "R.A. [arcmin]")
        img.axs[i][j].set_yticks([int(nx - k * nx / 5) for k in range(6)])
        img.axs[i][j].set_yticklabels(["%.1f" % (k * nx / 5. * pixel_scale) for k in range(6)])
        img.axs[i][j].set_xticks([int(k * ny / 5) for k in range(6)])
        img.axs[i][j].set_xticklabels(["%.1f" % (k * ny / 5. * pixel_scale) for k in range(6)])
img.save_img(parent_path + "shear_field.png")
img.show_img()

# the galactic parameters
h5f = h5py.File(parent_path + "param.hdf5", "w")

h5f["/FIELD/g"] = g[0]
h5f["/FIELD/g1"] = g[1]
h5f["/FIELD/g2"] = g[2]

for i in range(expo_num):

    ra = gal_coord[0,i*num_each_expo: (i+1)*num_each_expo]
    dec = gal_coord[1,i*num_each_expo: (i+1)*num_each_expo]
    print(ra.min(), ra.max(), ra.mean())
    print(dec.min(), dec.max(), dec.mean())

    h5f["/RA/expo_%d"%i] = ra
    h5f["/DEC/expo_%d"%i] = dec
    h5f["/RADIUS/expo_%d"%i] = gal_coord[2, i*num_each_expo: (i+1)*num_each_expo]

    g, g1, g2 = shear_field(ra, dec, shift, sigma, amplitude)
    h5f["/g/expo_%d"%i] = g
    h5f["/g1/expo_%d"%i] = g1
    h5f["/g2/expo_%d"%i] = g2
    print(g1.min(), g1.max())
    mag_i = tool_box.mag_generator(num_each_expo, 20, 23.5)
    flux_i = tool_box.mag_to_flux(mag_i)
    h5f["/MAG/expo_%d"%i] = mag_i
    h5f["/FLUX/expo_%d"%i] = flux_i
    print(flux_i.min(), flux_i.max())
h5f.close()





