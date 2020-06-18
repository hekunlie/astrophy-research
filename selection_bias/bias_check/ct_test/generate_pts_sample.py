from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
from astropy.io import fits
import h5py


psf_type = "Moffat"
psf_flux = 1
psf_scale = 4
stamp_size = 40
seed = 14432

noise_sig = 60

pts_num = 40
max_radius = 5
gal_flux = 25000

g1, g2 = 0.05, 0.05

total_num = 99

fq = Fourier_Quad(stamp_size, seed)


# psf_img = fq.cre_psf(psf_scale, psf_flux, psf_type)
# hdu = fits.PrimaryHDU(numpy.float32(psf_img))
# hdu.writeto("./imgs/psf.fits", overwrite=True)

# dt = numpy.linspace(-numpy.pi, numpy.pi, total_num)

pts = fq.ran_pts(pts_num, max_radius)
theta = numpy.random.uniform(0, 2*numpy.pi, 3600000)

h5f = h5py.File("./imgs/pts.hdf5","w")
h5f['/pts'] = numpy.float32(pts)
h5f['/theta'] = numpy.float32(theta)
h5f.close()
# for i in range(total_num):
#     pts_r = fq.rotate(pts,theta[i])
#     gal_img_nf = fq.convolve_psf(pts_r, psf_scale, gal_flux / pts_num, psf_type)
#     # img = Image_Plot()
#     # img.subplots(1,1)
#     # img.axs[0][0].imshow(gal_img_nf + fq.draw_noise(0,60))
#     # img.show_img()
#     # img.close_img()
#     hdu = fits.PrimaryHDU(numpy.float32(gal_img_nf))
#     hdu.writeto("./imgs/gal_%d.fits"%i,overwrite=True)