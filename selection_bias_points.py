from sys import path
path.append('/home/hklee/work/fourier_quad/')
import numpy
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import lsstetc
import time
import os
from multiprocessing import Pool
import pandas
import tool_box
# import galsim

def simu(paths_list, shear1, shear2, num_in_chip, magnitudes, proc_id):
    print('Process %d: simulation begin>>>>') % proc_id

    stamp_size = 56
    pixel_scale = 0.2
    psf_r = 4
    p_num = 45

    # the information from measurements
    cat_col = ["KSB_g1", "BJ_e1", "RG_e1", "FQ_G1", "fg1", "KSB_g2", "BJ_e2", "RG_e2", "FQ_G2", "fg2", "FG_N", "FQ_U",
               "FQ_V", 'KSB_R', 'BJ_R', 'RG_R', "area", "total_flux", "peak", "fnsr", "snr"]

    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
    noise_sig = prop.sigma_sky/20

    psf_in = Fourier_Quad().cre_psf(psf_r, stamp_size, 'Moffat')

    chips_num = len(paths_list)
    psf_pow = Fourier_Quad().pow_spec(psf_in)
    # psf_g = galsim.Image(psf_in)

    rim = Fourier_Quad().border(1, stamp_size)
    n = numpy.sum(rim)
    for path_tag in range(chips_num):
        t1 = time.time()
        chip_path = paths_list[path_tag]
        shear_tag, chip_name = chip_path.split('/')[3:5]
        chip_tag = chip_name.split('_')[2].split('.')[0]
        data_path = '/lmc/selection_bias/result/data/%s_gal_chip_%s.xlsx' % (shear_tag, chip_tag)

        # parameters
        g1_input = shear1[int(shear_tag)]
        g2_input = shear2[int(shear_tag)]
        tag = range(num_in_chip*int(chip_tag), num_in_chip*(1 + int(chip_tag)))
        mags = magnitudes[tag]
        data_matrix = numpy.zeros((num_in_chip, len(cat_col)))
        gal_pool = []
        for k in range(num_in_chip):
            gal_flux = prop.flux(mags[k])/p_num
            points = Fourier_Quad().ran_pos(num=p_num, radius=8, g=(g1_input, g2_input))[1]
            noise = numpy.random.normal(loc=0, scale=noise_sig, size=stamp_size ** 2).reshape(stamp_size, stamp_size)
            gal_final = Fourier_Quad().convolve_psf(pos=points, psf_scale=psf_r, imagesize=stamp_size,
                                                    flux=gal_flux, psf='Moffat')+noise

            gal_pool.append(gal_final)
            obj, flux, signalsq, peak = tool_box.stamp_detector(gal_final, noise_sig*1.5, stamp_size, stamp_size)
            snr = numpy.sqrt(signalsq) / noise_sig

            gpow = Fourier_Quad().pow_spec(gal_final)
            signal = numpy.sum(gpow[stamp_size/2-1:stamp_size/2+2, stamp_size/2-1:stamp_size/2+2])/9
            noise_level = numpy.sum(rim*gpow)/n
            fsnr = numpy.sqrt(signal/noise_level)


            # shear estimate

            # gal_g = galsim.Image(gal_final)
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

            noise_n = numpy.random.normal(loc=0., scale=noise_sig, size=stamp_size**2).reshape(stamp_size, stamp_size)
            mg1, mg2, mn, mu, mv = Fourier_Quad().shear_est(gal_final, psf_pow, stamp_size, noise_n, F=True)

            data_matrix[k, :] = 0, 0, 0, mg1, g1_input, 0, 0, 0, mg2, g2_input, mn, mu, mv, 0, 0, 0, \
                                len(obj), flux, peak, fsnr, snr

        data_df = pandas.DataFrame(data=data_matrix, columns=cat_col)
        data_df.to_excel(data_path)

        big_chip = Fourier_Quad().image_stack(gal_pool, stamp_size, 100)
        hdu = fits.PrimaryHDU(big_chip)
        hdu.writeto(chip_path, overwrite=True)
        t2 = time.time()
        print('Proc_%d: simulation %d/%d finished within %.2f<<<') % (proc_id, path_tag+1, chips_num, t2-t1)


if __name__ == '__main__':
    CPU_num = 16
    chip_num = 100
    total_num = 1000000
    num = total_num/chip_num

    data_files = os.listdir('/lmc/selection_bias/result/data/')
    if len(data_files) != 0:
        for name in data_files:
            os.remove('/lmc/selection_bias/result/data/' + name)
    for i in range(10):
        files_path = '/lmc/selection_bias/%d/' % i
        if not os.path.isdir(files_path):
            os.mkdir(files_path)
        else:
            files = os.listdir(files_path)
            if len(files) != 0:
                for name in files:
                    os.remove(files_path + name)

    chip_paths_pool = ['/lmc/selection_bias/%d/gal_chip_%s.fits' % (i, str(j).zfill(3)) for i in range(10)
                       for j in range(chip_num)]
    chip_paths_list = tool_box.task_distri(chip_paths_pool, CPU_num)

    # magnitude
    out_mags = numpy.load('/home/hklee/work/selection_bias/parameters/lsstmagsims.npz')['arr_0']

    # the input shear signal
    shear = numpy.load('/home/hklee/work/selection_bias/parameters/shear.npz')['arr_0']
    out_shear1 = shear[0]
    out_shear2 = shear[1]

    p = Pool()
    ts = time.time()
    for m in range(CPU_num):
        p.apply_async(simu, args=(chip_paths_list[m], out_shear1, out_shear2, num, out_mags, m,))
    p.close()
    p.join()
    # simu(chip_paths_list[0], out_shear1, out_shear2, num, out_mags, 0, 1,)
    te = time.time()
    print('Time consuming: %.2f') % (te - ts)
    os.system(" sudo python selection_bias_proce_dat.py 0 10000000 none 1")