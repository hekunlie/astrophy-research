import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import *
import matplotlib.pyplot as plt
import lsstetc
import time
from multiprocessing import Pool

def simulate(g1,g2,NO):

    print('Process %d: begin>>>>>>')%NO

    gal_num = 10000
    stamp_size = 60
    pixel_scale = 0.3

    mag_list = numpy.sort(numpy.loadtxt('lsstmagsims'))
    prop = lsstetc.ETC(band='i', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)

    psf  = galsim.Gaussian(half_light_radius=1.0)
    psf_img = psf.drawImage(nx=stamp_size, ny=stamp_size, scale=pixel_scale)

    ellip1 = numpy.random.normal(loc=0, scale=0.1, size=1000000)
    ellip2 = numpy.random.normal(loc=0, scale=0.1, size=1000000)

    for k in range(100):      # 100 chips
        kk = str(k).zfill(3)  # chip number
        ts = time.time()
        print('Process %d: Chip %s')%(NO,kk)

        chip_path = ''
        data_path = ''
        psf_path  = ''

        gal_pool   = []
        noise_pool = []

        snr_data   = numpy.zeros((gal_num,4))
        tag        = range(int(k * 10000), int((k + 1) * 10000))
        maglists   = mag_list[tag]
        ell1       = ellip1[tag]
        ell2       = ellip2[tag]

        for i in range(gal_num):
            e1 = ell1[i]
            e2 = ell2[i]
            mag= maglists[i]
            n = numpy.random.randint(1,3,1)[0]

            if n==1:
                r0  = numpy.random.randint(15,26)[0]/10.
                gal = galsim.Exponential(flux=1.0, half_light_radius=r0)

            elif n==2:
                r0    = numpy.random.randint(10, 13)[0]/10.
                bulge = galsim.DeVaucouleurs(flux=0.4, half_light_radius=r0)
                disk  = galsim.Exponential(flux=0.6,half_light_radius=2*r0)
                gal   = galsim.Add([bulge,disk])

            else:
                r0  = numpy.random.randint(15, 26)[0] / 10.
                gal = galsim.Gaussian(flux=1.0,half_light_radius=r0)


            gal = gal.shear(e1=e1,e2=e2)
            gal = gal.shear(g1=g1,g2=g2)
            gal = galsim.Convolve([psf,gal]) #the final galaxy profile

            gal_img, noise_img = prop.draw(gal,mag,i,noise=1)

            snr = prop.SNR(gal,mag)
            err = prop.err(gal,mag)

            snr_data[i,0:4] = snr, mag, prop.sigma_sky, err

            gal_pool.append(gal_img.array)
            noise_pool.append(noise_img.array)


        numpy.savetxt('/home/lmc/Downloads/snr.txt',snr_data)

        final = Fourier_Quad().image_stack(gal_pool,48,100)
        hdu = fits.PrimaryHDU(final)
        hdu.writeto('/home/lihekun/Downloads/gal_image.fits',overwrite=True)

        final = Fourier_Quad().image_stack(noise_pool,48,100)
        hdu = fits.PrimaryHDU(final)
        hdu.writeto('/home/lihekun/Downloads/noise_image.fits',overwrite=True)

        hdu = fits.PrimaryHDU(final)
        hdu.writeto('/home/lihekun/Downloads/psf_image.fits',overwrite=True)
        te = time.time()
        print('Process %d: Completed with time comsuming: %.2f')%(NO,te-ts)

if __name__=='__main__':
    shear1 = numpy.linspace(-0.005, 0.005, 11)
    shear2 = numpy.linspace(-0.005, 0.005, 11)
    numpy.random.shuffle(shear1)
    numpy.random.shuffle(shear2)
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