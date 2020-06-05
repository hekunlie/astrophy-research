from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
from astropy.io import fits
import h5py


psf_type = "Moffat"
psf_flux = 1
psf_scale = 4
stamp_size = 48
seed = 1412

noise_sig = 60

pts_num = 40
max_radius = 7
gal_flux = 6000

g1, g2 = 0, 0

nx,ny = 10,10

total_num = nx*ny

fq = Fourier_Quad(stamp_size, seed)

psf_img = fq.cre_psf(psf_scale, psf_flux, psf_type)
psf_pow = fq.pow_spec(psf_img)
fq.get_radius(psf_pow, 2)
print(fq.hlr)

shear_est_nf = numpy.zeros((total_num, 5))
shear_est_n = numpy.zeros((total_num, 5))
shear_est_ct = numpy.zeros((total_num, 5))
shear_est_ct_est = numpy.zeros((total_num, 5))
shear_est_nr = numpy.zeros((total_num, 5))

img_buffer = numpy.zeros((ny*stamp_size, nx*stamp_size))
img_buffer_n = numpy.zeros((ny*stamp_size, nx*stamp_size))

pts = fq.ran_pts(pts_num, max_radius)

dt = numpy.linspace(-numpy.pi, numpy.pi, total_num)
for i in range(total_num):

    noise = fq.draw_noise(0, noise_sig)
    noise_pow = fq.pow_spec(noise)

    noise_new = fq.draw_noise(0, noise_sig)
    noise_new_pow = fq.pow_spec(noise_new)

    noise_pow_ct = fq.pow_spec(noise + noise_new) - noise_new - noise_new_pow

    # noise residual power
    noise_pow_residual = fq.pow_spec(noise + noise_new) - noise_pow - noise_new_pow
    shear_est_nr[i] = fq.shear_est(noise_pow_residual, psf_pow, F=True)

    # noise free power
    pst_r = fq.rotate(pts, dt[i])
    gal_img_nf = fq.convolve_psf(pst_r, psf_scale, gal_flux/pts_num, psf_type)
    gal_pow_nf = fq.pow_spec(gal_img_nf)
    shear_est_nf[i] = fq.shear_est(gal_pow_nf, psf_pow, F=True)

    # noisy power
    gal_pow_n = fq.pow_spec(gal_img_nf + noise)
    gal_pow = gal_pow_n - noise_new_pow
    shear_est_n[i] = fq.shear_est(gal_pow, psf_pow, F=True)

    # cross term
    ct_pow = gal_pow - gal_pow_nf - noise_pow
    shear_est_ct[i] = fq.shear_est(ct_pow, psf_pow, F=True)

    # est cross term
    gal_img_new = gal_img_nf + noise + noise_new
    gal_img_new_pow = fq.pow_spec(gal_img_new)
    ct_pow_est = gal_img_new_pow - gal_pow_n - noise_new_pow - noise_pow_ct
    shear_est_ct_est[i] = fq.shear_est(ct_pow_est, psf_pow, F=True)

    # stack
    iy, ix = divmod(i, nx)
    img_buffer[iy*stamp_size:(iy+1)*stamp_size, ix*stamp_size:(ix+1)*stamp_size] = gal_img_nf
    img_buffer_n[iy*stamp_size:(iy+1)*stamp_size, ix*stamp_size:(ix+1)*stamp_size] = gal_img_nf+noise


h5f = h5py.File("./sample_data.hdf5", "w")
h5f["/noise_free"] = shear_est_nf
h5f["/noise_residual"] = shear_est_nr
h5f["/noisy"] = shear_est_n
h5f["/cross_term"] = shear_est_ct
h5f["/cross_term_est"] = shear_est_ct_est
h5f.close()

print(shear_est_ct_est.shape)
hdu = fits.PrimaryHDU(img_buffer)
hdu.writeto("./sample.fits",overwrite=True)
hdu = fits.PrimaryHDU(img_buffer_n)
hdu.writeto("./sample_noisy.fits",overwrite=True)