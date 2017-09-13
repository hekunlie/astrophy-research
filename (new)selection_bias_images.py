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

def simu(gal_paths_list, shear1_in, shear2_in, num_in_chip, e1_in, e2_in, proc_id):
    print('Process %d: begin>>>>') % proc_id

    stamp_size = 80
    pixel_scale = 0.2
    col = ['morphology', 'mag', 'snr', 'noise_sigma']

    prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)

    psf_o = galsim.Gaussian(half_light_radius=0.8, flux=1.0)
    psf = psf_o.shear(e1=0.081, e2=-0.066)


    for i in range(len(gal_paths_list)):
        chip_path = gal_paths_list[i]
        shear_tag, chip_name = chip_path.split('/')[3:5]
        info_path = '/lmc/selection_bias/%s/gal_info_%s.xlsx' % (shear_tag, chip_name.split('_')[2].split('.')[0])
        g1_input = shear1_in[int(shear_tag)]
        g2_input = shear2_in[int(shear_tag)]
        pool = []
        for k in range(num_in_chip):

            morpho = numpy.random.randint(1, 4, 1)[0]
            ra = radius[i]
            if morpho == 1:
                gal = galsim.Exponential(flux=gal_flux, half_light_radius=ra)

            elif morpho == 2:
                rb = gal_radius_sb[rs_tag]
                rd = gal_radius_sd[rs_tag]
                rs_tag += 1
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
            snr_data[i, 0:4] = morpho, mag, snr, prop.sigma_sky
            gal_pool.append(gal_img)

    pass


if __name__=='__main__':
    CPU_num = 20
    chip_num = 250
    total_num = 1000000
    stamp_size = 80

    for i in range(10):
        files_path = '/lmc/selection_bias/%d/'%i
        if not os.path.isdir(files_path):
            os.mkdir(files_path)

    chip_paths_pool = ['/lmc/selection_bias/%d/gal_chip_%s.fits'%(i, str(j).zfill(2)) for i in range(10) for j in range(chip_num)]
    chip_paths_list = tool_box.task_distri(chip_paths_pool, CPU_num)

    arr = numpy.load('/lmc/selection_bias/shear.npz')
    shear1 = arr['arr_0']
    shear2 = arr['arr_1']
    mags = numpy.load('/home/hklee/work/selection_bias/parameters/lsstmagsims.npz')['arr_0']
    gal_rad = numpy.load('/home/hklee/work/selection_bias/parameters/gal_radius.npz')['arr_0']
    sersic_rad = numpy.load('/home/hklee/work/selection_bias/parameters/sersic_rd.npz')
    sersic_bulge = sersic_rad['arr_0']
    sersic_disk = sersic_rad['arr_1']
    e1e2 = numpy.load('/home/hklee/work/selection_bias/parameters/e1e2.npz')
    ie1 = e1e2['arr_0']
    ie2 = e1e2['arr_1']

    p = Pool()
    t1 = time.time()
    for m in range(CPU_num):
        p.apply_async(simu, args=(chip_paths_list,shear1, shear2, ie1, ie2, m,))
    p.close()
    p.join()
    # simu(shear1[0], shear2[0], ie1, ie2, gal_ra, sersic_rd, mags, 0)
    t2 = time.time()
    print('Time consuming: %.2f') % (t2 - t1)
    os.system('python selection_bias_est.py')