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
import tool_box

def simu(paths_list, psf_in, shear1_in, shear2_in, num_in_chip, e1_in, e2_in, magni, gal_radius, sersic_rb,
         sersic_rd, proc_id, est_switch):
    print('Process %d: simulation begin>>>>') % proc_id

    stamp_size = 80
    pixel_scale = 0.2

    # the initial information
    info_col = ['morphology', 'mag', 'snr', 'hla_flux', 'peak_val', 'peak_y', 'peak_x']
    # the information from measurements
    cat_col = ["KSB_g1", "BJ_e1", "RG_e1", "FQ_G1", "fg1", "KSB_g2", "BJ_e2", "RG_e2", "FQ_G2", "fg2", "FG_N", "FQ_U",
               "FQ_V", 'KSB_R', 'BJ_R', 'RG_R', "SNR_ORI", "hl_flux", "peak"]

    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
    noise_sig = prop.sigma_sky

    psf_o = galsim.Gaussian(half_light_radius=0.8, flux=1.0)
    psf = psf_o.shear(e1=0.081, e2=-0.066)

    chips_num = len(paths_list)
    psf_pow = Fourier_Quad().pow_spec(psf_in)
    psf_g = galsim.Image(psf_in)

    for path_tag in range(chips_num):
        t1 = time.time()
        chip_path = paths_list[path_tag]
        shear_tag, chip_name = chip_path.split('/')[3:5]
        chip_tag = chip_name.split('_')[2].split('.')[0]
        info_path = '/lmc/selection_bias/%s/gal_info_%s.xlsx' % (shear_tag, chip_tag)
        data_path = '/lmc/selection_bias/result/data/%s_gal_chip_%s.xlsx' % (shear_tag, chip_tag)

        # parameters
        g1_input = shear1_in[int(shear_tag)]
        g2_input = shear2_in[int(shear_tag)]
        tag = range(num_in_chip*int(chip_tag), num_in_chip*(1 + int(chip_tag)))
        mags = magni[tag]
        ellip1 = e1_in[tag]
        ellip2 = e2_in[tag]
        radius = gal_radius[tag]
        gal_radius_sb = sersic_rb[tag]
        gal_radius_sd = sersic_rd[tag]

        gal_pool = []
        snr_data = numpy.zeros((num_in_chip, len(info_col)))
        data_matrix = numpy.zeros((num_in_chip, len(cat_col)))
        for k in range(num_in_chip):
            e1 = ellip1[k]
            e2 = ellip2[k]
            gal_flux = prop.flux(mags[k])
            morpho = numpy.random.randint(1, 4, 1)[0]
            ra = radius[k]
            if morpho == 1:
                gal = galsim.Exponential(flux=gal_flux, half_light_radius=ra)

            elif morpho == 3:
                rb = gal_radius_sb[k]
                rd = gal_radius_sd[k]
                bulge = galsim.Sersic(half_light_radius=rb, n=3.5)
                disk = galsim.Sersic(half_light_radius=rd, n=1.5)
                gal = bulge * 0.3 + disk * 0.7
                gal = gal.withFlux(gal_flux)

            else:
                gal = galsim.Gaussian(flux=gal_flux, half_light_radius=ra)

            gal_s = gal.shear(e1=e1, e2=e2)
            gal_g = gal_s.shear(g1=g1_input, g2=g2_input)
            gal_c = galsim.Convolve([psf, gal_g])
            gal_img, snr = prop.draw(gal_c, add_noise=1)
            gal_pool.append(gal_img)
            hl_flux, maximum, cords = Fourier_Quad().get_radius_new(gal_img, 2., stamp_size)[2:5]
            snr_data[k] = morpho, mags[k], snr, hl_flux, maximum, cords[0], cords[1]

            # shear estimate
            if est_switch == 1:
                gal_g = galsim.Image(gal_img)

                res_k = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='KSB', strict=False)
                ksb_g1 = res_k.corrected_g1
                ksb_g2 = res_k.corrected_g2
                ksb_r = res_k.resolution_factor

                res_b = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='BJ', strict=False)
                bj_e1 = res_b.corrected_e1
                bj_e2 = res_b.corrected_e2
                bj_r = res_b.resolution_factor

                res_r = galsim.hsm.EstimateShear(gal_g, psf_g, shear_est='REGAUSS', strict=False)
                re_e1 = res_r.corrected_e1
                re_e2 = res_r.corrected_e2
                re_r = res_r.resolution_factor

                noise = numpy.random.normal(loc=0., scale=noise_sig, size=stamp_size**2).reshape(stamp_size, stamp_size)
                G1, G2, N, U, V = Fourier_Quad().shear_est(gal_img, psf_pow, stamp_size, noise, F=True, N=True)

                data_matrix[k, :] = ksb_g1, bj_e1, re_e1, G1, g1_input, ksb_g2, bj_e2, re_e2, G2, g2_input, \
                                    N, U, V, ksb_r, bj_r, re_r, snr, hl_flux, maximum

        info_df = pandas.DataFrame(data=snr_data, columns=info_col)
        info_df.to_excel(info_path)

        if est_switch == 1:
            data_df = pandas.DataFrame(data=data_matrix, columns=cat_col)
            data_df.to_excel(data_path)

        big_chip = Fourier_Quad().image_stack(gal_pool, stamp_size, 50)
        hdu = fits.PrimaryHDU(big_chip)
        hdu.writeto(chip_path, overwrite=True)
        t2 = time.time()
        print('Proc_%d: simulation %d/%d finished within %.2f<<<') % (proc_id, path_tag+1, chips_num, t2-t1)


if __name__ == '__main__':
    CPU_num = 16
    chip_num = 250
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

    # the input shear signal
    shear = numpy.load('/lmc/selection_bias/shear.npz')
    out_shear1 = shear['arr_0']
    out_shear2 = shear['arr_1']

    # magnitude
    out_mags = numpy.load('/home/hklee/work/selection_bias/parameters/lsstmagsims.npz')['arr_0']

    # galactic radius
    out_gal_rad = numpy.load('/home/hklee/work/selection_bias/parameters/gal_radius.npz')['arr_0']

    # the input Sersic parameters
    out_sersic_rad = numpy.load('/home/hklee/work/selection_bias/parameters/sersic_rd.npz')
    out_sersic_bulge = out_sersic_rad['arr_0']
    out_sersic_disk = out_sersic_rad['arr_1']

    # the input ellipticity
    out_e1e2 = numpy.load('/home/hklee/work/selection_bias/parameters/e1e2.npz')
    out_ie1 = out_e1e2['arr_0']
    out_ie2 = out_e1e2['arr_1']

    # the psf
    if not os.path.exists('/lmc/selection_bias/psf.fits'):
        out_prop = lsstetc.ETC(band='r', pixel_scale=0.2, stamp_size=80, nvisits=180)
        out_psf_flux = out_prop.flux(20.5)
        out_psf = galsim.Gaussian(half_light_radius=0.8, flux=1.0)
        out_psf = out_psf.shear(e1=0.081, e2=-0.066).withFlux(out_psf_flux)
        psf_m = out_prop.draw(out_psf, add_noise=1)[0]
    else:
        psf_m = fits.open('/lmc/selection_bias/psf.fits')[0].data

    p = Pool()
    ts = time.time()
    for m in range(CPU_num):
        p.apply_async(simu, args=(chip_paths_list[m], psf_m, out_shear1, out_shear2, num, out_ie1, out_ie2, out_mags,
                                  out_gal_rad, out_sersic_bulge, out_sersic_disk, m, 1,))
    p.close()
    p.join()
    # simu(chip_paths_list[0], psf_m, out_shear1, out_shear2, num, out_ie1, out_ie2, out_mags, out_gal_rad, out_sersic_bulge, out_sersic_disk, 0, 1,)
    te = time.time()
    print('Time consuming: %.2f') % (te - ts)