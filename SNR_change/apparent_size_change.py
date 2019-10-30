import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
from subprocess import Popen
import numpy
from plot_tool import Image_Plot
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import tool_box
import galsim


size = int(argv[1])
psf_r = float(argv[2])
e = float(argv[3])
btr = float(argv[4])
ra = float(argv[5])
seed = int(argv[6])
file_tag = int(argv[7])

pts_source = 0

num = 11
pixel_scale = 0.187

flux = [tool_box.mag_to_flux(21.5)]

noise_sig = 60

detect_thresh = 2

fq = Fourier_Quad(size, seed)
fq_p = Fourier_Quad(size, 17060)
# all the images are added by the same noise
noise = fq.draw_noise(0, noise_sig)

psf = galsim.Moffat(beta=3.5, fwhm=psf_r, flux=1.0, trunc=psf_r*3)

shear_beta = numpy.linspace(0, numpy.pi, num)
input_g = numpy.linspace(-0.06, 0.06, num)

pts_num = 100
rand_pts = fq_p.ran_pts(num=pts_num, radius=10, ellip=0.8)

gal_pool = []

if pts_source == 0:
    # rng = galsim.BaseDeviate(12300000)
    # knot = galsim.randwalk.RandomWalk(npoints=60, half_light_radius=ra, flux=1, rng=rng)
    bulge = galsim.Sersic(half_light_radius=ra, n=4, trunc=4.5 * ra)  # be careful
    disk = galsim.Sersic(scale_radius=ra, n=1, trunc=4.5 * ra)  # be careful
    gal = bulge * btr + disk * (1 - btr)  # + ktr*knot
    gal = gal.shear(e1=e, e2=0)  # beta=0.*galsim.degrees)
    gal_f = gal.withFlux(flux[0])
    #
    # rng = galsim.BaseDeviate(12300000)
    # gal = galsim.randwalk.RandomWalk(npoints=200, half_light_radius=ra, flux=flux[k], rng=rng)
    # gal_f = gal.shear(e1=e, e2=0)

    gal_c = galsim.Convolve([gal_f, psf])
    img_0 = galsim.ImageD(size, size)
    gal_c.drawImage(image=img_0, scale=pixel_scale)
    ori_gal_img = img_0.array + noise
# random walk galaxy
else:
    img_0 = fq.convolve_psf(rand_pts, 4, flux[0] / pts_num, "Moffat")
    ori_gal_img = img_0 + noise

gal_pool.append(ori_gal_img)

# the sheared galaxy
for i in range(num):
    if pts_source == 0:
        gal_s = gal_f.shear(g1=input_g[i], g2=0)  # beta=shear_beta[i]*galsim.radians)
        gal_s_c = galsim.Convolve([gal_s, psf])
        img_s = galsim.ImageD(size, size)
        gal_s_c.drawImage(image=img_s, scale=pixel_scale)
        shear_gal_img = img_s.array + noise
    else:
        rand_pts_s = fq.shear(rand_pts, g1=input_g[i], g2=0)
        img_s = fq.convolve_psf(rand_pts_s, 4, flux[k] / pts_num, "Moffat")
        shear_gal_img = img_s + noise

    gal_pool.append(shear_gal_img)
    path = img_path + "gal_shear_%d_%d.fits" % (k, i)
    hdu = fits.PrimaryHDU(shear_gal_img)
    hdu.writeto(path, overwrite=True)

