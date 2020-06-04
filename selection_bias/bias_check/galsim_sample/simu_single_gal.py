import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import time
from mpi4py import MPI
import h5py
import tool_box


parent_path = argv[1]

seed_step = 1
e1 = 0
e2 = 0.6
sersic_idx = 0.31
scale_radius = 0.3
gal_flux = 8000
noise_sig = 60
pixel_scale = 0.187
stamp_size = 48


psf = galsim.Moffat(beta=3.5, fwhm=0.7, flux=1.0, trunc=1.4).shear(e1=0.1, e2=0.)

psf_img = galsim.ImageD(stamp_size, stamp_size)
psf.drawImage(image=psf_img, scale=pixel_scale)
psf_arr = numpy.float32(psf_img.array)
hdu = fits.PrimaryHDU(psf_arr)
psf_path = parent_path + '/psf.fits'
hdu.writeto(psf_path, overwrite=True)


gal = galsim.Sersic(scale_radius=scale_radius, n=sersic_idx, trunc=4.5 * scale_radius, flux=1.0)


g1_input = [0, -0.04,      0,   0.04,  -0.02,        0, 0.02,      0,  0.02]
g2_input = [0,        0, 0.04,  -0.04,       0,  -0.02,      0, 0.02, -0.02]


for shear_id in range(len(g1_input)):

    g1 = g1_input[shear_id]
    g2 = g2_input[shear_id]

    chip_path = parent_path + "/gal_%d.fits" %shear_id

    gal_e = gal.shear(e1=e1, e2=e2).withFlux(gal_flux)
    gal_s = gal_e.shear(g1=g1, g2=g2)
    gal_c = galsim.Convolve([gal_s, psf])
    gal_img = galsim.ImageD(stamp_size, stamp_size)
    gal_c.drawImage(image=gal_img, scale=pixel_scale)
    hdu = fits.PrimaryHDU(numpy.float32(gal_img.array))
    hdu.writeto(chip_path, overwrite=True)


