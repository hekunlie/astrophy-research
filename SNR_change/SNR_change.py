import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
from astropy.io import fits
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
import galsim
import time

size = int(argv[1])
psf_r = float(argv[2])
e = float(argv[3])
btr = float(argv[4])
ra = float(argv[5])


seed = numpy.random.randint(0,100000,1)[0]
num = 11
pixel_scale = 0.187

flux = numpy.array([tool_box.mag_to_flux(21.5),tool_box.mag_to_flux(22.2), tool_box.mag_to_flux(24.1)])
sig = 60
print(sig)

fq = Fourier_Quad(size, seed)
# all the images are added by the same noise
noise = fq.draw_noise(0, sig)

psf = galsim.Moffat(beta=3.5, fwhm=psf_r, flux=1.0, trunc=psf_r*3)

shear_beta = numpy.linspace(0, numpy.pi, num)
input_g = numpy.linspace(-0.06, 0.06, num)
img_path = os.getcwd() + "/imgs/"
for k in range(len(flux)):

    pool = []
    rng = galsim.BaseDeviate(12300000)
    #knot = galsim.randwalk.RandomWalk(npoints=60, half_light_radius=ra, flux=1, rng=rng)
    # bulge = galsim.Sersic(half_light_radius=ra, n=4, trunc=4.5 * ra)  # be careful
    # disk = galsim.Sersic(scale_radius=ra, n=1, trunc=4.5 * ra)  # be careful
    # gal = bulge * btr + disk * (1 - btr) #+ ktr*knot
    # gal = gal.shear(e1=e, e2=0)#beta=0.*galsim.degrees)
    # gal_f = gal.withFlux(flux[k])

    # rng = galsim.BaseDeviate(12300000)
    # gal = galsim.randwalk.RandomWalk(npoints=200, half_light_radius=ra, flux=flux[k], rng=rng)
    # gal_f = gal.shear(e1=e, e2=0)


    # gal_c = galsim.Convolve([gal_f, psf])
    # img_0 = galsim.ImageD(size, size)
    # gal_c.drawImage(image=img_0, scale=pixel_scale)
    # gal_img = img_0.array + noise

    rand_pts = fq.ran_pos(100, 9)
    img_0 = fq.convolve_psf(rand_pts, 4, flux[k]/100, "Moffat")

    gal_img = img_0 + noise
    path = img_path + "gal0_%d.fits" %k
    hdu = fits.PrimaryHDU(gal_img)
    hdu.writeto(path, overwrite=True)

    for i in range(num):

        # gal_s = gal_f.shear(g1=input_g[i], g2=0)#beta=shear_beta[i]*galsim.radians)
        # gal_s_c = galsim.Convolve([gal_s, psf])
        # img_s = galsim.ImageD(size, size)
        # gal_s_c.drawImage(image=img_s, scale=pixel_scale)
        # gal_s_img = img_s.array + noise

        rand_pts_s = fq.shear(rand_pts, g1=input_g[i], g2=0)
        img_s = fq.convolve_psf(rand_pts_s, 4, flux[k]/100, "Moffat")
        gal_s_img = img_s + noise

        path = img_path + "gal_%d_%d.fits"%(k,i)
        hdu = fits.PrimaryHDU(gal_s_img)
        hdu.writeto(path, overwrite=True)

time.sleep(10)
title = "S_%d\\PR_%.2f\\e1_%.2f\\BTR_%.2f\\GR_%.2f\\seed_%d"%(size, psf_r, e, btr, ra, seed)
cmd = "python SNR_change_plot.py %d %d %s %d"%(size, num, title, len(flux))
os.system(cmd)
