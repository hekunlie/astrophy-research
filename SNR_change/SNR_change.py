import matplotlib
matplotlib.use("Agg")
from sys import path
path.append('/home/hkli/work/fourier_quad/')
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from Fourier_Quad import Fourier_Quad
import tool_box
import lsstetc
from subprocess import Popen
import os
from sys import argv
import galsim
import copy

size = int(argv[1])
psf_r = float(argv[2])
ra = float(argv[3])
e = float(argv[4])
btr = float(argv[5])
seed = int(argv[6])


num = 17
pixel_scale = 0.2
markers = ['p', 'x', 's', '8']
colors = ['red', 'orange', 'green', 'deepskyblue', 'b', 'k']

prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=size, nvisits=180)
flux = numpy.array([prop.flux(22.5), prop.flux(23.8),  prop.flux(24.9)])
sig = prop.sigma_sky
print(sig)

fq = Fourier_Quad(size, seed)
# all the images are added by the same noise
noise = fq.draw_noise(0, sig)

psf = galsim.Moffat(beta=3.5, scale_radius=psf_r, flux=1.0)

shear_beta = numpy.linspace(0, numpy.pi, num)
input_g = numpy.linspace(-0.06, 0.06, num)
for k in range(len(flux)):

    pool = []
    bulge = galsim.Sersic(half_light_radius=0.6 * ra, n=4, trunc=4.5 * ra)  # be careful
    disk = galsim.Sersic(half_light_radius=ra, n=1, trunc=4.5 * ra)  # be careful
    gal = bulge * btr + disk * (1 - btr)
    gal = gal.shear(e1=e, e2=0)#beta=0.*galsim.degrees)
    gal_f = gal.withFlux(flux[k])

    gal_c = galsim.Convolve([gal_f, psf])
    img_0 = galsim.ImageD(size, size)
    gal_c.drawImage(image=img_0, scale=pixel_scale)
    gal_img = img_0.array + noise
    path = "/home/hkli/work/sex/imgs/gal0_%d.fits" %k
    hdu = fits.PrimaryHDU(gal_img)
    hdu.writeto(path, overwrite=True)

    for i in range(num):

        gal_s = gal_f.shear(g1=input_g[i], g2=0)#beta=shear_beta[i]*galsim.radians)
        gal_s_c = galsim.Convolve([gal_s, psf])
        img_s = galsim.ImageD(size, size)
        gal_s_c.drawImage(image=img_s, scale=pixel_scale)
        gal_s_img = img_s.array + noise

        path = "/home/hkli/work/sex/imgs/gal_%d_%d.fits"%(k,i)
        hdu = fits.PrimaryHDU(gal_s_img)
        hdu.writeto(path, overwrite=True)
title = "xy_%d\\pr_%.2f\\gr_%.2f\\e1%.2f\\btr_%.2f\\seed_%d"%(size, psf_r,ra, e, btr, seed)
cmd = "python SNR_change_plot.py %d %d %s"%(size, num, title)
os.system(cmd)
