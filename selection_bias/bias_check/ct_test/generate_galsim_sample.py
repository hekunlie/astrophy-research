from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
from astropy.io import fits
import h5py
import galsim
import os

sersic_idx = 0.4
scale_radius = 0.4
gal_flux = 8000
noise_sig = 60


pixel_scale = 0.187
stamp_size = 48

psf = galsim.Moffat(beta=3.5, fwhm=0.7, flux=1.0, trunc=1.4)
psf_img = galsim.ImageD(stamp_size, stamp_size)
psf.drawImage(image=psf_img, scale=pixel_scale)
psf_arr = numpy.float32(psf_img.array)
hdu = fits.PrimaryHDU(psf_arr)
psf_path = "./imgs/psf.fits"
hdu.writeto(psf_path, overwrite=True)

gal = galsim.Sersic(scale_radius=scale_radius, n=sersic_idx, trunc=4.5 * scale_radius, flux=1.0)

# nx = 9
# ny = 5
#
# big_img = numpy.zeros((ny*stamp_size, nx*stamp_size),dtype=numpy.float32)
# e1 = numpy.linspace(-0.7, 0.7, nx)
# e2 = numpy.linspace(-0.6, 0.6, ny)
#
# for i in range(ny):
#     for j in range(nx):
#
#         tag = i*nx + j
#
#         data_path = "/mnt/perc/hklee/bias_check/ct_test/scatter_test/cpsf_galsim/imgs/%d"%tag
#         if not os.path.exists(data_path):
#             os.makedirs(data_path)
#
#         gal_e = gal.shear(e1=e1[j], e2=e2[i]).withFlux(gal_flux)
#         gal_c = galsim.Convolve([gal_e, psf])
#         gal_img = galsim.ImageD(stamp_size, stamp_size)
#         gal_c.drawImage(image=gal_img, scale=pixel_scale)
#
#         gal_img = numpy.float32(gal_img.array)
#         big_img[i*stamp_size:(i+1)*stamp_size, j*stamp_size:(j+1)*stamp_size] = gal_img
#
#         hdu = fits.PrimaryHDU(gal_img)
#         hdu.writeto("./imgs/gal_%d.fits"%tag,overwrite=True)
#
# hdu = fits.PrimaryHDU(big_img)
# hdu.writeto("./imgs/gal_all.fits",overwrite=True)
# big_img = big_img + numpy.random.normal(0,noise_sig,(ny*stamp_size, nx*stamp_size))
# hdu = fits.PrimaryHDU(big_img)
# hdu.writeto("./imgs/gal_all_noisy.fits",overwrite=True)

nx = 9
ny = 5

big_img = numpy.zeros((ny*stamp_size, nx*stamp_size),dtype=numpy.float32)
e1 = numpy.linspace(0.1,0.7,ny)
theta = numpy.linspace(0,90,nx)
for i in range(ny):
    for j in range(nx):

        tag = i*nx + j

        data_path = "/mnt/perc/hklee/bias_check/ct_test/scatter_test/cpsf_galsim/imgs/%d"%tag
        if not os.path.exists(data_path):
            os.makedirs(data_path)

        gal_e = gal.shear(e1=e1[i], e2=0).rotate(theta[j]*galsim.degrees).withFlux(gal_flux)
        gal_c = galsim.Convolve([gal_e, psf])
        gal_img = galsim.ImageD(stamp_size, stamp_size)
        gal_c.drawImage(image=gal_img, scale=pixel_scale)

        gal_img = numpy.float32(gal_img.array)
        big_img[i*stamp_size:(i+1)*stamp_size, j*stamp_size:(j+1)*stamp_size] = gal_img

        hdu = fits.PrimaryHDU(gal_img)
        hdu.writeto("./imgs/gal_%d.fits"%tag,overwrite=True)

hdu = fits.PrimaryHDU(big_img)
hdu.writeto("./imgs/gal_all.fits",overwrite=True)
big_img = big_img + numpy.random.normal(0,noise_sig,(ny*stamp_size, nx*stamp_size))
hdu = fits.PrimaryHDU(big_img)
hdu.writeto("./imgs/gal_all_noisy.fits",overwrite=True)
