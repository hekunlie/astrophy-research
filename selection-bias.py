from sys import path
path.append('/home/hklee/codes/')
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import lsstetc
import time
import os
from multiprocessing import Pool
import pandas

def simulate(g1, g2, ellip1, ellip2, gal_radius, gal_radius_s, noise_imgs, mag_list, process_id):

    print('Process %d: begin>>>>')%process_id
    gal_num = 10000
    chip_num = 100
    stamp_size = 80
    pixel_scale = 0.2
    col = ['morphology', 'mag', 'snr', 'noise_sigma']
    label = range(0, gal_num)
    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)

    psf_o = galsim.Gaussian(half_light_radius=1.0)
    psf = psf_o.shear(e1=0.08, e2=-0.06)
    psf_noise = numpy.random.normal(loc=0, scale=380.8645, size=stamp_size**2).reshape(stamp_size, stamp_size)
    psf_img = prop.draw(psf, 21, psf_noise, add_noise=1)[0]

    ahead = '/lmc/selection_bias/%d/' %process_id
    if not os.path.isdir(ahead):
        os.mkdir(ahead)
    for k in range(chip_num):
        kk = str(k).zfill(2)
        ts = time.time()
        print('Process %d: Simulate Chip %s')%(process_id, kk)

        gal_chip_path = ahead + 'gal_chip_%s.fits'%kk
        data_path = ahead + 'gal_info_%s.xlsx'%kk
        psf_path = ahead + 'psf.fits'

        gal_pool = []
        snr_data = numpy.zeros((gal_num, 4))
        tag = range(int(k * gal_num), int((k + 1) * gal_num))
        mag_piece = mag_list[tag]
        ell1 = ellip1[tag]
        ell2 = ellip2[tag]
        radius = gal_radius[tag]
        noise_stamps = noise_imgs[tag]
        rs_tag = 0
        for i in range(gal_num):
            e1 = ell1[i]
            e2 = ell2[i]
            mag = mag_piece[i]
            morpho = numpy.random.randint(1, 4, 1)[0]
            ra = radius[i]
            rs = gal_radius_s[rs_tag]

            if morpho==1:
                gal = galsim.Exponential(flux=1.0, half_light_radius=ra)

            elif morpho==2:
                rs_tag += 1
                rb = ra/2
                rd = ra/2 + rs
                bulge = galsim.Sersic(half_light_radius=rb, n=3.5)
                disk = galsim.Sersic(half_light_radius=rd, n=1.5)
                gal = bulge*0.3+disk*0.7

            else:
                gal = galsim.Gaussian(flux=1.0, half_light_radius=ra)

            gal_s = gal.shear(e1=e1, e2=e2)
            gal_g = gal_s.shear(g1=g1, g2=g2)
            gal_c = galsim.Convolve([psf, gal_g])
            noise_stamp = noise_stamps[i]
            gal_img, snr = prop.draw(gal_c, mag, noise_stamp, add_noise=1)
            snr_data[i, 0:4] = morpho, mag, snr, prop.sigma_sky
            gal_pool.append(gal_img)

        df = pandas.DataFrame(data=snr_data, index=label, columns=col)
        df.to_excel(data_path)

        gal_chip = Fourier_Quad().image_stack(gal_pool, stamp_size, 100)
        hdu = fits.PrimaryHDU(gal_chip)
        hdu.writeto(gal_chip_path, overwrite=True)

        hdu = fits.PrimaryHDU(psf_img)
        hdu.writeto(psf_path, overwrite=True)

        te = time.time()
        print('Process %d: Simulation completed with time consuming: %.2f'%(process_id, te-ts))

if __name__=='__main__':
    arr = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = arr['arr_0']
    shear2 = arr['arr_1']
    mags = numpy.load('/home/hklee/work/selection_bias/parameters/lsstmagsims.npz')['arr_0']
    noise = numpy.load('/home/hklee/work/selection_bias/parameters/lsstmagsims.npz')['arr_0']
    gal_ra = numpy.load('/home/hklee/work/selection_bias/parameters/gal_radius.npz')['arr_0']
    sersic_rd = numpy.load('/home/hklee/work/selection_bias/parameters/sersic_rd.npz')['arr_0']
    e1e2 = numpy.load('/home/hklee/work/selection_bias/parameters/e1e2.npz')
    ie1 = e1e2['arr_0']
    ie2 = e1e2['arr_1']
    numpy.random.shuffle(mags)
    p = Pool()
    t1 = time.time()
    for m in range(len(shear1)):
       ig1 = shear1[m]
       ig2 = shear2[m]
       p.apply_async(simulate, args=(ig1, ig2, ie1, ie2, gal_ra, sersic_rd, noise, mags, m,))
    p.close()
    p.join()
    t2 = time.time()
    print('Time consuming: %.2f') % (t2 - t1)
    os.system('python selection_bias_est.py')

