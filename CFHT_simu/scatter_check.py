import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import time
from mpi4py import MPI
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

source = argv[1]

envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['selection_bias', "%s_path" % source, '1'], ['selection_bias', "%s_path_result" % source, '1'],
                ['selection_bias', "%s_path_para" % source, '1'], ['selection_bias', "%s_path_log" % source, '1']]
path_items = tool_box.config(envs_path, ['get', 'get', 'get', 'get'], get_contents)
total_path, result_path, para_path, log_path = path_items


stamp_size = 48
stamp_col = 100
noise_sig = 60
pixel_scale = 0.187
stamp_num = 10000

mag = 21
ra = 0.8
btr = 0.4
e1 = 0.6
e2 = -0.7
gal_flux = tool_box.mag_to_flux(mag)

psf = galsim.Moffat(beta=3.5, fwhm=0.7, flux=1.0, trunc=2)

gal_pool = []

rng = numpy.random.RandomState(rank+123*rank+1)
fq = Fourier_Quad(stamp_size, 123)
ny,nx = stamp_col*stamp_size, stamp_col*stamp_size

bulge = galsim.Sersic(half_light_radius=ra, n=4, trunc=4.5 * ra, flux=1.0)  # be careful
disk = galsim.Sersic(scale_radius=ra, n=1, trunc=4.5 * ra, flux=1.0)  # be careful
gal = bulge * btr + disk * (1 - btr)
gal_e = gal.shear(e1=e1, e2=e2)
gal_f = gal_e.withFlux(gal_flux)
gal_c = galsim.Convolve([gal_f, psf])
img = galsim.ImageD(stamp_size, stamp_size)
gal_c.drawImage(image=img, scale=pixel_scale)
gal_arr = img.array

for i in range(stamp_num):
    gal_pool.append(gal_arr)

noise_img = rng.normal(0, noise_sig, nx * ny).reshape((ny, nx))
big_chip = fq.stack(gal_pool, stamp_col) + noise_img
big_chip = numpy.float32(big_chip)
hdu = fits.PrimaryHDU(big_chip)
chip_path = total_path + "scatter/%d.fits"%rank
hdu.writeto(chip_path, overwrite=True)
te = time.clock()
print(te-ts)