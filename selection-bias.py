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

    print('Process %d: Simulation begin>>>>>>')%NO
    t11 = time.time()
    gal_num = 2000
    chip_num = 1
    total_num = chip_num*gal_num*1000
    stamp_size = 80
    pixel_scale = 0.2
    col = ['mophology', 'mag', 'snr', 'noise_sigma']
    label = range(0, gal_num)
    mag_list = numpy.loadtxt('lsstmagsims')
    numpy.random.shuffle(mag_list)
    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)

    psf  = galsim.Gaussian(half_light_radius=0.8)
    #psf  = psf.shear(e1=0.05,e2=-0.03)
    psf_img = psf.drawImage(nx=stamp_size, ny=stamp_size, scale=pixel_scale).array

    ellip = numpy.random.normal(loc=0, scale=0.18, size=total_num)
    idx1 = ellip < 0.8
    idx2 = ellip > 0.
    num = len(ellip[idx1 & idx2])
    ellip = ellip[idx1 & idx2]
    dnum = int(num - total_num / 2)
    if dnum != 0:
        if dnum < 0:
            e_add = numpy.random.uniform(0, 0.8, -dnum)
            ellip = numpy.append(ellip, e_add)
        else:
            ellip = ellip[0:int(total_num / 2)]
    ellip = numpy.append(-ellip, ellip)
    numpy.random.shuffle(ellip)
    while True:
        theta = 4*numpy.random.uniform(0,1.00000001,total_num)*numpy.pi
        ellip1 = ellip*numpy.cos(theta)
        ellip2 = ellip*numpy.sin(theta)
        if numpy.abs(numpy.mean(ellip2))<1.e-5 and numpy.abs(numpy.mean(ellip1))<1.e-5:
            break

    ahead = '/lmc/selection_bias/%d/' %NO
    if not os.path.isdir(ahead):
        os.mkdir(ahead)
    seed_ori = numpy.random.randint(0, total_num*10, total_num)
    t22 =time.time()
    print('the first part %.2f'%(t22-t11))

    for k in range(chip_num):
        kk = str(k).zfill(2)
        ts = time.time()
        print('Process %d: Simulate Chip %s')%(NO,kk)

        gal_chip_path = ahead + 'gal_chip_%s.fits'%kk
        noise_chip_path = ahead + 'noise_chip_%s.fits'%kk
        data_path = ahead + 'gal_info_%s.xlsx'%kk
        psf_path  = ahead + 'psf.fits'

        gal_pool   = []
        noise_pool = []

        snr_data  = numpy.zeros((gal_num,4))
        tag       = range(int(k * gal_num), int((k + 1) * gal_num))
        mag_piece = mag_list[tag]
        ell1      = ellip1[tag]
        ell2      = ellip2[tag]
        r0 = numpy.random.uniform(0.4,1.2,gal_num)
        seed_k = seed_ori[tag]
        tm = time.time()
        tgau = 0
        ng = 0
        tser = 0
        ns =0
        texp = 0
        ne=0
        for i in range(gal_num):
            tgals =time.time()
            e1 = ell1[i]
            e2 = ell2[i]
            mag= mag_piece[i]
            mopho = numpy.random.randint(1,4,1)[0]
            r = r0[i]

            if mopho==1:
                gal = galsim.Exponential(flux=1.0, half_light_radius=r)

            elif mopho==2:
                rb    = r0[i]/2.5 + 0.3
                rd    = r0[i]/2.0 + 1.1
                bulge = galsim.Sersic(half_light_radius=rb,n=3.5)
                disk  = galsim.Sersic(scale_radius=rd,n=2.0)
                gal   = bulge*0.3+disk*0.7

            else:
                gal = galsim.Gaussian(flux=1.0,half_light_radius=r)

            gal_s = gal.shear(e1=e1,e2=e2)
            gal_g = gal_s.shear(g1=g1,g2=g2)
            gal_c = galsim.Convolve([psf,gal_g]) #the final galaxy profile
            seed = seed_k[i]
            gal_img, noise_img, snr = prop.draw(gal_c,mag,seed,add_noise=1)
            #gal_img = prop.draw(gal_c, mag, seed, add_noise=None)
            snr_data[i,0:4] = mopho, mag, snr, prop.sigma_sky
            gal_pool.append(gal_img.array)
            noise_pool.append(noise_img.array)
            tgale = time.time()
            delt = tgale-tgals
            if mopho ==3:
                tgau+=delt
                ng+=1
            elif mopho==2:
                tser+=delt
                ns+=1
            else:
                texp+=delt
                ne+=1
        print('gaussian %.2f, sersic %.2f, exponential %.2f'%(tgau/ng,tser/ns,texp/ne))
        df = pandas.DataFrame(data= snr_data, index = label, columns=col)
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
        print('Process %d: Simulation completed with time comsuming: %.2f,%.2f'%(NO,te-tm,tm-ts))

if __name__=='__main__':
    arr = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = arr['arr_0']
    shear2 = arr['arr_1']
    p = Pool()
    t1 = time.time()
    # for m in range(len(shear1)):
    #    g1 = shear1[m]
    #    g2 = shear2[m]
    #    p.apply_async(simulate, args=(g1,g2,m,))
    # p.close()
    # p.join()
    simulate(shear1[1],shear2[1],1)
    t2 = time.time()
    print('Time comsuming: %.2f') % (t2 - t1)
    #os.system('python selection_bias_est.py')

