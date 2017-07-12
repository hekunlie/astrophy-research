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
    col = ['mophology', 'snr', 'mag', 'noise_sigma', 'err']
    label = range(0, gal_num)
    mag_list = numpy.sort(numpy.loadtxt('lsstmagsims'))
    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)

    psf  = galsim.Gaussian(half_light_radius=1.2)
    psf  = psf.shear(e1=0.05,e2=-0.03)

    psf_img = psf.drawImage(nx=stamp_size, ny=stamp_size, scale=pixel_scale).array
    ellip1_gener = numpy.abs(numpy.random.normal(loc=0, scale=0.15, size=500000))
    idx = ellip1_gener >0.6
    ellip1_gener[idx] = 0.6
    ellip1 = numpy.append(-ellip1_gener,ellip1_gener)

    ellip2_gener = numpy.abs(numpy.random.normal(loc=0, scale=0.15, size=500000))
    idx = ellip2_gener >0.6
    ellip2_gener[idx] = 0.6
    ellip2 = numpy.append(-ellip2_gener,ellip2_gener)

    numpy.random.shuffle(ellip1)
    numpy.random.shuffle(ellip2)

    ahead = '/lmc/selection_bias/%d/' %NO
    if not os.path.isdir(ahead):
        os.mkdir(ahead)

    for k in range(100):      # 100 chips
        kk = str(k).zfill(2)      # chip number
        ts = time.time()
        print('Process %d: Chip %s')%(NO,kk)

        gal_chip_path = ahead + 'gal_chip_%s.fits'%kk
        noise_chip_path = ahead + 'noise_chip_%s.fits'%kk
        data_path = ahead + 'gal_info_%s.xlsx'%kk
        psf_path    = ahead + 'psf.fits'

        gal_pool   = []
        noise_pool = []

        snr_data    = numpy.zeros((gal_num,5))
        tag             = range(int(k * gal_num), int((k + 1) * gal_num))
        mag_piece = mag_list[tag]
        ell1       = ellip1[tag]
        ell2       = ellip2[tag]
        r0 = numpy.random.random(gal_num)
        idx = r0 <0.3
        r0[idx]=0.3
        for i in range(gal_num):
            e1 = ell1[i]
            e2 = ell2[i]
            mag= mag_piece[i]
            mopho = numpy.random.randint(1,4,1)[0]
            r = r0[i] * 1.2 + 0.3

            if mopho==1:
                gal = galsim.Exponential(flux=1.0, half_light_radius=r)

            elif mopho==2:
                rd     = r0[i]+0.5
                bulge = galsim.Sersic(half_light_radius=r,n=3.5)
                disk  = galsim.Sersic(scale_radius=rd,n=1.5)
                gal   = bulge*0.3+disk*0.7

            else:
                gal = galsim.Gaussian(flux=1.0,half_light_radius=r)

            gal_s = gal.shear(e1=e1,e2=e2)
            gal_g = gal_s.shear(g1=g1,g2=g2)
            gal_c = galsim.Convolve([psf,gal_g]) #the final galaxy profile

            gal_img, noise_img, snr = prop.draw(gal_c,mag,i,add_noise=1)
            snr_data[i,0:5] = mopho, snr, mag, prop.sigma_sky, 2.5/numpy.log(10)/snr
            gal_pool.append(gal_img.array)
            noise_pool.append(noise_img.array)


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
    # shear1 = numpy.linspace(-0.005, 0.005, 11)
    # shear2 = numpy.linspace(-0.005, 0.005, 11)
    # numpy.random.shuffle(shear1)
    arr = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = arr['arr_0']
    shear2 = arr['arr_1']
    p = Pool()
    t1 = time.time()
    for m in range(1,11):
       g1 = shear1[m]
       g2 = shear2[m]
       p.apply_async(simulate, args=(g1,g2,m,))
    p.close()
    p.join()

    t2 = time.time()
    print('Time comsuming: %.2f')%(t2-t1)
