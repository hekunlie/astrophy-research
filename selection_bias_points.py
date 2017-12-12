from __future__ import division
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import numpy
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import lsstetc
import time
import tool_box
from mpi4py import MPI
# import galsim
import logging


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

t_s = time.clock()

time.sleep(rank*0.1)

with open("/home/hklee/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "total_data" in path:
        total_path = path.split("=")[1]
    elif "result" in path:
        result_path = path.split("=")[1]
    elif "parameters" in path:
        para_path = path.split("=")[1]
    elif "logs_path" in path:
        logs_path = path.split("=")[1]

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logfile = logs_path + '%d_log.txt' %rank

lf = logging.FileHandler(logfile, 'w')
form = logging.Formatter('%(asctime)s - %(message)s')
lf.setFormatter(form)
logger.addHandler(lf)

stamp_size = 54
pixel_scale = 0.2
psf_r = 4
psf_model = "Moffat"
p_num = 45
chip_num = 500
total_num = 5000000
num_in_chip = int(total_num/chip_num)

fq = Fourier_Quad(stamp_size, rank*12344+31121)

# distribute jobs
chip_paths_pool = [total_path + '%d/gal_chip_%s.fits' % (i, str(j).zfill(3)) for i in range(10)
                   for j in range(chip_num)]
chip_paths_list = tool_box.task_distri(chip_paths_pool, cpus)[rank]

# LSST
prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
noise_sig = prop.sigma_sky

# magnitude
magnitude = numpy.load(para_path + 'lsstmagsims.npz')['arr_0']

# the input shear signal
shear = numpy.load(para_path + 'shear.npz')['arr_0']
shear1 = shear[0]
shear2 = shear[1]

# the information from measurements
cat_col = ["KSB_g1", "BJ_e1", "RG_e1", "FQ_G1", "fg1", "KSB_g2", "BJ_e2", "RG_e2", "FQ_G2", "fg2", "FG_N", "FQ_U",
           "FQ_V", 'KSB_R', 'BJ_R', 'RG_R', "area", "total_flux", "peak", "fnsr", "snr", "ori_snr", "flag"]

psf_in = fq.cre_psf(psf_r, psf_model)
psf_pow = fq.pow_spec(psf_in)

chips_num_indiv = len(chip_paths_list)
rim = fq.border(1)
n = numpy.sum(rim)

if rank == 0:
    seqs = ["stamp size:\t%d\n"%stamp_size, "pixel scale:\t%.1f\n"%pixel_scale, "chips number:\t%d\n"%chip_num,
            "total galaxy number:\t%d\n"%total_num, "PSF model:\t%s\n"%psf_model, "PSF scale:\t%d\n"%psf_r,
            "Noise sigma:\t%.2f\n"%noise_sig, "Magnitudes:\t%.2f ~ %.2f\n"%(numpy.min(magnitude), numpy.max(magnitude))]

    with open(result_path + "data/paras.dat", "w") as para:
        para.writelines(seqs)

logger.info("%d's process: gets %d (%d) chips" % (rank, chips_num_indiv, len(chip_paths_pool)))

for path_tag in range(chips_num_indiv):

    logger.info("%d's process: %d's (%d) chips starts..." % (rank, path_tag+1, chips_num_indiv))

    t1 = time.time()
    chip_path = chip_paths_list[path_tag]
    shear_tag, chip_name = chip_path.split('/')[3:5]
    chip_tag = chip_name.split('_')[2].split('.')[0]
    data_path = result_path + 'data/%s_gal_chip_%s.npz' % (shear_tag, chip_tag)

    # parameters
    g1_input = shear1[int(shear_tag)]
    g2_input = shear2[int(shear_tag)]
    tag = range(num_in_chip*int(chip_tag), num_in_chip*(1 + int(chip_tag)))
    mags = magnitude[tag]
    data_matrix = numpy.zeros((num_in_chip, len(cat_col)))
    gal_pool = []

    logger.info("%d's process: %d's (%d) chips prepared the parameters" % (rank, path_tag + 1, chips_num_indiv))

    for k in range(num_in_chip):
        gal_flux = prop.flux(mags[k])/p_num
        points = fq.ran_pos(num=p_num, radius=9, g=(g1_input, g2_input))[1]
        noise = fq.noise(0, noise_sig)
        gal_final = fq.convolve_psf(pos=points, psf_scale=psf_r, flux=gal_flux, psf=psf_model) + noise

        gal_pool.append(gal_final)
        obj, flux, signalsq, peak, flag = tool_box.stamp_detector(gal_final, noise_sig*2, stamp_size, stamp_size)
        snr = numpy.sqrt(signalsq) / noise_sig
        ori_snr = flux/numpy.sqrt(len(obj))/noise_sig

        gpow = fq.pow_spec(gal_final)
        if numpy.max(gpow) == gpow[int(stamp_size/2), int(stamp_size/2)]:
            signal = gpow[int(stamp_size/2), int(stamp_size/2)]
        else:
            signal = numpy.sum(gpow[int(stamp_size/2-1):int(stamp_size/2+2), int(stamp_size/2-1):int(stamp_size/2+2)])/9
        noise_level = numpy.sum(rim*gpow)/n
        fsnr = numpy.sqrt(signal/noise_level)

        # shear estimate

        noise_n = fq.noise(0, noise_sig)
        mg1, mg2, mn, mu, mv = fq.shear_est(gal_final, psf_pow, noise_n, F=True)

        data_matrix[k, :] = 0, 0, 0, mg1, g1_input, 0, 0, 0, mg2, g2_input, mn, mu, mv, 0, 0, 0, \
                            len(obj), flux, peak, fsnr, snr, ori_snr, flag

    logger.info("%d's process: %d's (%d) chips's writing to files." % (rank, path_tag+1, chips_num_indiv))

    # data_df = pandas.DataFrame(data=data_matrix, columns=cat_col)
    # data_df.to_excel(data_path)
    numpy.savez(data_path, data_matrix)
    if int(chip_tag) < 30:
        big_chip = fq.image_stack(gal_pool, 100)
        hdu = fits.PrimaryHDU(big_chip)
        hdu.writeto(chip_path, overwrite=True)
    t2 = time.time()

    logger.info("%d's process: %d's (%d) chips finished within %.2f." % (rank, path_tag+1, chips_num_indiv, t2-t1))

t_e = time.clock()
logger.info("%d's process: finished within %.2f." % (rank, t_e-t_s))
