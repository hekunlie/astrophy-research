import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import lsstetc
import time
import os
from multiprocessing import Pool
import pandas

def simulate(g1, g2, NO):

    print('Process %d: begin>>>>>>')%NO

    gal_num = 10000
    stamp_size = 60
    pixel_scale = 0.2

    mag_list = numpy.sort(numpy.loadtxt('lsstmagsims'))
    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)

    psf  = galsim.Gaussian(half_light_radius=0.7)
    psf  = psf.shear(e1=0.01,e2=-0.03)
    psf_img = psf.drawImage(nx=stamp_size, ny=stamp_size, scale=pixel_scale).array

    ellip1 = numpy.random.normal(loc=0, scale=0.1, size=1000000)
    ellip2 = numpy.random.normal(loc=0, scale=0.1, size=1000000)

    ahead = '/lmc/selection_bias/%d/' %NO
    if not os.path.isdir(ahead):
        os.mkdir(ahead)

    for k in range(100):      # 100 chips
        kk = str(k).zfill(2)      # chip number
        ts = time.time()
        print('Process %d: Chip %s')%(NO,kk)

        gal_chip_path = ahead + 'gal_chip_%s.fits'%kk
        noise_chip_path = ahead + 'nosie_chip_%s.fits'%kk
        data_path = ahead + 'gal_info_%s.xlsx'%kk
        psf_path    = ahead + 'psf.fits'

        gal_pool   = []
        noise_pool = []

        snr_data    = numpy.zeros((gal_num,5))
        tag             = range(int(k * 10000), int((k + 1) * 10000))
        mag_piece = mag_list[tag]
        ell1       = ellip1[tag]
        ell2       = ellip2[tag]

        for i in range(gal_num):
            e1 = ell1[i]
            e2 = ell2[i]
            mag= mag_piece[i]
            mopho = int(numpy.random.randint(1,4,1)[0])

            if mopho==1:
                r0  = numpy.random.randint(15,27,1)[0]/10.
                gal = galsim.Exponential(flux=1.0, half_light_radius=r0)

            elif mopho==2:
                r0    = numpy.random.randint(10, 14,1)[0]/10.
                bulge = galsim.DeVaucouleurs(flux=0.4, half_light_radius=r0)
                disk  = galsim.Exponential(flux=0.6,half_light_radius=2*r0)
                gal   = galsim.Add([bulge,disk])

            else:
                r0  = numpy.random.randint(15, 27,1)[0] / 10.
                gal = galsim.Gaussian(flux=1.0,half_light_radius=r0)

            if e1**2+e2**2>1:
                e1 = e1/2
                e2 = e2/2
            gal = gal.shear(e1=e1,e2=e2)
            gal = gal.shear(g1=g1,g2=g2)
            gal = galsim.Convolve([psf,gal]) #the final galaxy profile

            gal_img, noise_img = prop.draw(gal,mag,i,noise=1)

            snr = prop.SNR(gal,mag)
            err = prop.err(gal,mag)

            snr_data[i,0:5] = mopho, snr, mag, prop.sigma_sky, err

            gal_pool.append(gal_img.array)
            noise_pool.append(noise_img.array)

        label = range(0,10000)
        col = ['snr','mag','noise_sigma','err']
        df = pandas.DataFrame(data= snr_data,index = label,columns=col)
        df.to_excel(data_path)

        gal_chip = Fourier_Quad().image_stack(gal_pool,stamp_size,100)
        hdu = fits.PrimaryHDU(gal_chip)
        hdu.writeto(gal_chip_path,overwrite=True)

        noise_chip = Fourier_Quad().image_stack(noise_pool,stamp_size,100)
        hdu = fits.PrimaryHDU(noise_chip)
        hdu.writeto(noise_chip_path,overwrite=True)

        hdu = fits.PrimaryHDU(psf_img)
        hdu.writeto(psf_path,overwrite=True)

        te = time.time()
        print('Process %d: Completed with time comsuming: %.2f')%(NO,te-ts)

if __name__=='__main__':
    shear1 = numpy.linspace(-0.005, 0.005, 11)
    shear2 = numpy.linspace(-0.005, 0.005, 11)
    numpy.random.shuffle(shear1)
    numpy.savez('/lmc/selection_bias/shear',shear1,shear2)
    p = Pool()
    t1 = time.time()
    for m in range(11):
        g1 = shear1[m]
        g2 = shear2[m]
        p.apply_async(simulate, args=(g1,g2,m,))
    p.close()
    p.join()
    t2 = time.time()
    print('Time comsuming: %.2f')%(t2-t1)