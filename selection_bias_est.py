from sys import path
path.append('/home/hklee/codes/')
from Fourier_Quad import  Fourier_Quad
import numpy
from astropy.io import fits
import os
from multiprocessing import Pool
import pandas
import galsim
import time

def est_shear(m, g1, g2):
    stamp_size = 80
    chip_num = 100
    ahead = '/lmc/selection_bias/%s/'%m
    respath = '/lmc/selection_bias/result/data/'
    if not os.path.isdir(respath):
        os.makedirs(respath)
    psf = fits.open(ahead+'psf.fits')[0].data
    col =  ["KSB_g1", "BJ_e1", "RG_e1", "FQ_G1", "FG_N", "fg1", "KSB_g2", "BJ_e2", "RG_e2", "FQ_G2", "FG_N", "fg2", "FQ_U", "FQ_V", "SNR_ORI"]
    for k in range(chip_num):
        ts = time.time()
        kk = str(k).zfill(2)
        print('Process %d: Shear estimation chip %s beginning>>>>')%(m, kk)
        gal_path = ahead + 'gal_chip_%s.fits'%kk
        gal_img  = fits.open(gal_path)[0].data
        # noise_path = ahead + 'noise_chip_%s.fits'%kk
        # noise_img = fits.open(noise_path)[0].data
        cat_path  = ahead + 'gal_info_%s.xlsx'%kk
        cat_data  = pandas.read_excel(cat_path).values[:, 2]

        gal_pool = Fourier_Quad().divide_stamps(gal_img, stamp_size)
        # noise_pool = Fourier_Quad().divide_stamps(noise_img, stamp_size)
        gal_index = []
        data_matrix = numpy.zeros((len(gal_pool), len(col)))
        for j in range(len(gal_pool)):
            gg = str(j).zfill(4)
            idx = kk+'_%s'%gg
            gal_index.append(idx)
            gal = gal_pool[j]
            noise = numpy.random.normal(loc=0., scale=380.4, size=stamp_size**2).reshape(stamp_size, stamp_size)
            gal_g = galsim.Image(gal)
            psf_g = galsim.Image(psf)
            res_k = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='KSB', strict=False)
            ksb_g1 = res_k.corrected_g1
            ksb_g2 = res_k.corrected_g2
            res_b = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='BJ', strict=False)
            bj_e1 = res_b.corrected_e1
            bj_e2 = res_b.corrected_e2
            res_r = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='REGAUSS', strict=False)
            re_e1 = res_r.corrected_e1
            re_e2 = res_r.corrected_e2
            G1, G2, N, U, V = Fourier_Quad().shear_est(gal, psf, stamp_size, noise, N=True)
            data_matrix[j, :] = ksb_g1, bj_e1, re_e1, G1, N, g1, ksb_g2, bj_e2, re_e2, G2, N, g2, U, V, cat_data[j]

            #data_matrix[j,:] = 0, 0, 0, G1, N, g1, 0, 0, 0, G2, N, g2, U, V,cat_data[j]
        df = pandas.DataFrame(data_matrix, index=gal_index, columns=col)
        df.columns.name = 'Chip&NO'
        res_path = respath+'%d_chip_%s.xlsx'%(m, kk)
        df.to_excel(res_path)
        te =time.time()
        print('Process %d: shear estimation chip %s complete with time consuming %.2f')%(m, kk, te-ts)

if __name__=='__main__':
    shear = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = shear['arr_0']
    shear2 = shear['arr_1']
    p = Pool()
    t1 = time.time()
    for m in range(len(shear1)):
       g1 = shear1[m]
       g2 = shear2[m]
       p.apply_async(est_shear, args=(m, g1, g2,))
    p.close()
    p.join()
    t2 = time.time()
    #est_shear(0,shear1[0],shear2[0])
    print('Time comsuming: %.2f') % (t2 - t1)
    #os.system('python selection_bias_proce_dat.py')
