import numpy
from sys import path, argv
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
path.append('%s/work/mylib/'%my_home)
import plot_tool
import matplotlib.pyplot as plt
from astropy.io import fits
import h5py


area_id = int(argv[1])
tag = int(argv[2])
radius_e = float(argv[3])

scale = 0.15

h5f = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_grid.hdf5","r")
shape = h5f["/background/w_%d"%area_id].attrs["grid_shape"]
ny, nx = shape[0], shape[1]
ra = h5f["/foreground/w_%d/RA"%area_id].value[tag]
dec = h5f["/foreground/w_%d/DEC"%area_id].value[tag]

ra_bin = h5f["/background/w_%d/RA_bin"%area_id].value
dec_bin = h5f["/background/w_%d/DEC_bin"%area_id].value

ra_min, ra_max = ra_bin.min(), ra_bin.max()
dec_min, dec_max = dec_bin.min(), dec_bin.max()

h5f.close()
print(ra,dec,ra_bin[1] - ra_bin[0]-scale, dec_bin[1] - dec_bin[0]-scale)


mask = fits.open("mask.fits")[0].data
target_blocks = []
for i in range(ny):
    for j in range(nx):
        if mask[i,j] > -1:
            m,n = divmod(mask[i,j], nx)
            target_blocks.append((m,n))
            print(m,n,dec_bin[m], dec_bin[m+1], ra_bin[n], ra_bin[n+1])
print(len(target_blocks))
img = plot_tool.Image_Plot(fig_x=20, fig_y=int(20.*ny/nx))
img.plot_img(1,1)

for i in range(ny + 1):
    img.axs[0][0].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linewidth=0.5)
    for j in range(nx + 1):
        img.axs[0][0].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black", linewidth=0.5)

for blks in target_blocks:
    blk_x, blks_y = ra_bin[blks[1]] + 0.5 * scale, dec_bin[blks[0]] + 0.5 * scale
    img.axs[0][0].scatter(blk_x, blks_y, s=5, c="blue")
img.axs[0][0].scatter(ra, dec, s=10, c="red")

theta = numpy.linspace(0,2*numpy.pi,10000)
re_y, re_x = dec + radius_e * numpy.sin(theta), ra + radius_e * numpy.cos(theta)
ey1 = re_y >= dec_min
ey2 = re_y <= dec_max
ex1 = re_x >= ra_min
ex2 = re_x <= ra_max
idx2 = ey1 & ey2 & ex1 & ex2

img.axs[0][0].scatter(re_x[idx2], re_y[idx2], c="coral", s=0.03)
limit = None
if limit:
    x_min, x_max = max(ra - 1.5*radius_e,ra_min), min(ra + 1.5*radius_e, ra_max)
    y_min, y_max = max(dec - 1.5*radius_e, dec_min), min(dec + 1.5*radius_e, dec_max)

    img.axs[0][0].set_xlim(x_min, x_max)
    img.axs[0][0].set_ylim(y_min, y_max)
img.save_img("radius.png")
