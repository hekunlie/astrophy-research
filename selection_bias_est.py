from sys import path
path.append('/home/hklee/work/fourier_quad/')
from Fourier_Quad import Fourier_Quad
import numpy
from astropy.io import fits
import os
from multiprocessing import Pool
import pandas
import galsim
import time
import lsstetc
import tool_box

def shear_est(chip_list, psf_in, shear1_in, shear2_in, noise_sig, size, proc_id):
    print('Proc_%d: begin>>>')%proc_id

    # col = ["KSB_g1", "BJ_e1", "RG_e1", "FQ_G1", "fg1", "KSB_g2", "BJ_e2", "RG_e2", "FQ_G2", "fg2",
    #        "FG_N", "FQ_U", "FQ_V",  'KSB_R', 'BJ_R', 'RG_R', "SNR_ORI"]

    psf_pow = Fourier_Quad().pow_spec(psf_in)
    # psf_g = galsim.Image(psf_in)

    total_chips = len(chip_list)
    for i in range(total_chips):
        chip_path = chip_list[i]
        shear_tag, chip_name = chip_path.split('/')[3:5]
        # info_path = '/lmc/selection_bias/%s/gal_info_%s.xlsx' %(shear_tag, chip_name.split('_')[2].split('.')[0])

        # g1_input = shear1_in[int(shear_tag)]
        # g2_input = shear2_in[int(shear_tag)]

        gals = fits.open(chip_path)[0].data
        gal_pool = Fourier_Quad().divide_stamps(gals, size)
        # snr = pandas.read_excel(info_path).values[:, 2]

        # data_matrix = numpy.zeros((len(gal_pool), len(col)))
        data_path = '/lmc/selection_bias/result/data/' + shear_tag + '_' + chip_name.split('.')[0] + '.xlsx'
        ts = time.time()
        data1 = numpy.zeros((len(gal_pool), 1))
        data2 = numpy.zeros((len(gal_pool), 2))
        for k in range(len(gal_pool)):
            gal = gal_pool[k]
            # gal_g = galsim.Image(gal)
            #
            # res_k = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='KSB', strict=False)
            # ksb_g1 = res_k.corrected_g1
            # ksb_g2 = res_k.corrected_g2
            # ksb_r = res_k.resolution_factor
            #
            # res_b = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='BJ', strict=False)
            # bj_e1 = res_b.corrected_e1
            # bj_e2 = res_b.corrected_e2
            # bj_r = res_b.resolution_factor
            #
            # res_r = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='REGAUSS', strict=False)
            # re_e1 = res_r.corrected_e1
            # re_e2 = res_r.corrected_e2
            # re_r = res_r.resolution_factor

            noise = numpy.random.normal(loc=0., scale=noise_sig, size=size**2).reshape(size, size)
            mg1, mg2, mn, mu, mv = Fourier_Quad().shear_est(gal, psf_pow, size, noise, F=True)

            # data_matrix[k, :] = ksb_g1, bj_e1, re_e1, mg1, g1_input, ksb_g2, bj_e2, re_e2, mg2, g2_input, mn, mu, mv, ksb_r, bj_r, re_r, snr[k]
            data1[k] = mg1
            data2[k] = mg2
        df = pandas.read_excel(data_path)
        df["FQ_G1"] = data1
        df["FQ_G2"] = data2
        df.to_excel(data_path)
        te = time.time()
        print('Proc_%d: (%d/%d) complete within time %.2f s') % (proc_id, i+1, total_chips, te-ts)

if __name__=='__main__':
    CPU_num = 16
    chip_num = 200
    total_num = 1000000
    pixel_scale = 0.2
    stamp_size = 58

    result_path = '/lmc/selection_bias/result/data/'
    if not os.path.isdir(result_path):
        os.makedirs(result_path)

    shear = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = shear['arr_0']
    shear2 = shear['arr_1']

    chip_paths_pool = ['/lmc/selection_bias/%d/gal_chip_%s.fits' % (i, str(j).zfill(3)) for i in range(10)
                       for j in range(chip_num)]
    chip_paths_list = tool_box.task_distri(chip_paths_pool, CPU_num)

    # psf = fits.open('/lmc/selection_bias/psf.fits')[0].data
    psf = Fourier_Quad().cre_psf(4, stamp_size, 'Moffat')
    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
    noise_sigma = prop.sigma_sky/8

    p = Pool()
    t1 = time.time()
    for m in range(CPU_num):
       p.apply_async(shear_est, args=(chip_paths_list[m], psf, shear1, shear2, noise_sigma, stamp_size, m))
    p.close()
    p.join()
    t2 = time.time()
    # shear_est(chip_paths_list[0], psf, shear1, shear2, noise_sigma, stamp_size, 0)
    print('Time comsuming: %.2f') % (t2 - t1)