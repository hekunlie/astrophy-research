from sys import path
path.append('/home/hklee/work/fourier_quad/')
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import lsstetc
import time
import os
from multiprocessing import Pool
import pandas

def simulate(g1, g2, ellip1, ellip2, gal_radius, gal_radius_s, mag_list, process_id):

    print('Process %d: begin>>>>')%process_id
    gal_num = 10000
    chip_num = 100
    stamp_size = 80
    pixel_scale = 0.2
    col = ['morphology', 'mag', 'snr', 'noise_sigma']
    label = range(0, gal_num)
    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)

    psf_o = galsim.Gaussian(half_light_radius=0.8, flux=1.0)
    psf = psf_o.shear(e1=0.081, e2=-0.066)
    psf_noise = numpy.random.normal(loc=0, scale=prop.sigma_sky, size=stamp_size**2).reshape(stamp_size, stamp_size)
    psf_img = psf.withFlux(prop.flux(20.5)).drawImage(nx=stamp_size, ny=stamp_size, scale=pixel_scale).array
    psf_img = psf_img + psf_noise

    ahead = '/lmc/selection_bias/%d/' %process_id
    if not os.path.isdir(ahead):
        os.mkdir(ahead)

    psf_path = ahead + 'psf.fits'
    hdu = fits.PrimaryHDU(psf_img)
    hdu.writeto(psf_path, overwrite=True)

    sersic_rb = gal_radius_s['arr_0']
    sersic_rd = gal_radius_s['arr_1']
    rs_tag = 0
    for k in range(chip_num):
        kk = str(k).zfill(2)
        ts = time.time()
        print('Process %d: Simulate Chip %s')%(process_id, kk)

        gal_chip_path = ahead + 'gal_chip_%s.fits'%kk
        data_path = ahead + 'gal_info_%s.xlsx'%kk

        gal_pool = []
        snr_data = numpy.zeros((gal_num, 4))
        tag = range(int(k * gal_num), int((k + 1) * gal_num))
        mag_piece = mag_list[tag]
        ell1 = ellip1[tag]
        ell2 = ellip2[tag]
        radius = gal_radius[tag]
        for i in range(gal_num):
            e1 = ell1[i]
            e2 = ell2[i]
            mag = mag_piece[i]
            gal_flux = prop.flux(mag)
            morpho = numpy.random.randint(1, 4, 1)[0]
            ra = radius[i]
            if morpho == 1:
                gal = galsim.Exponential(flux=gal_flux, half_light_radius=ra)

            elif morpho == 2:
                rb = sersic_rb[tag]
                rd = sersic_rd[tag]
                rs_tag += 1
                bulge = galsim.Sersic(half_light_radius=rb, n=3.5)
                disk = galsim.Sersic(half_light_radius=rd, n=1.5)
                gal = bulge*0.3 + disk*0.7
                gal = gal.withFlux(gal_flux)

            else:
                gal = galsim.Gaussian(flux=gal_flux, half_light_radius=ra)

            gal_s = gal.shear(e1=e1, e2=e2)
            gal_g = gal_s.shear(g1=g1, g2=g2)
            gal_c = galsim.Convolve([psf, gal_g])
            gal_img, snr = prop.draw(gal_c, add_noise=1)
            snr_data[i, 0:4] = morpho, mag, snr, prop.sigma_sky
            gal_pool.append(gal_img)

        df = pandas.DataFrame(data=snr_data, index=label, columns=col)
        df.to_excel(data_path)

        gal_chip = Fourier_Quad().image_stack(gal_pool, stamp_size, 100)
        hdu = fits.PrimaryHDU(gal_chip)
        hdu.writeto(gal_chip_path, overwrite=True)

        te = time.time()
        print('Process %d: Simulation completed with time consuming: %.2f'%(process_id, te-ts))

if __name__=='__main__':
    arr = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = arr['arr_0']
    shear2 = arr['arr_1']
    mags = numpy.load('/home/hklee/work/selection_bias/parameters/lsstmagsims.npz')['arr_0']
    gal_rad = numpy.load('/home/hklee/work/selection_bias/parameters/gal_radius.npz')['arr_0']
    sersic_rad = numpy.load('/home/hklee/work/selection_bias/parameters/sersic_rd.npz')
    e1e2 = numpy.load('/home/hklee/work/selection_bias/parameters/e1e2.npz')
    ie1 = e1e2['arr_0']
    ie2 = e1e2['arr_1']
    p = Pool()
    t1 = time.time()
    for m in range(len(shear1)):
       ig1 = shear1[m]
       ig2 = shear2[m]
       p.apply_async(simulate, args=(ig1, ig2, ie1, ie2, gal_rad, sersic_rad, mags, m,))
    p.close()
    p.join()
    #simulate(shear1[0], shear2[0], ie1, ie2, gal_ra, sersic_rd, mags, 0)
    t2 = time.time()
    print('Time consuming: %.2f') % (t2 - t1)
    os.system('python selection_bias_est.py')

