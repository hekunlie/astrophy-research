import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import lsstetc
import time
from mpi4py import MPI
import h5py
import logging
import tool_box

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

ini_path = "%s/work/envs/envs.dat"%my_home
path_items = tool_box.congif(ini_path,['get'], [['selection_bias', "dimmerm_path_para", '1']])

total_path = "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer_m/PDF/"
result_path = "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer_m/PDF/result/"
para_path = path_items[0]
pic_path = "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer_m/PDF/pic/"
log_path = "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer_m/PDF/log/"


logger = logging.getLogger()
logger.setLevel(logging.INFO)
logfile = log_path + '%d_log.dat' %rank

lf = logging.FileHandler(logfile, 'w')
form = logging.Formatter('%(asctime)s - %(message)s')
lf.setFormatter(form)
logger.addHandler(lf)

stamp_size = 90
pixel_scale = 0.2
shear_num = 1
chips_num = int(500/int(cpus/shear_num))
seed = rank*344 + 1121
chip_s_id, shear_id = divmod(rank, shear_num)

fq = Fourier_Quad(stamp_size, seed)

shear_cata = para_path + "shear.npz"
shear = numpy.load(shear_cata)
g1 = 0.06
g2 = -0.06

paras = para_path + "para_%d.hdf5"%shear_id
f = h5py.File(paras,'r')
e1s = f["/e1"].value
e2s = f["/e2"].value
radius = f["/radius"].value
flux = f["/flux"].value
fbt = f['/btr'].value
f.close()

prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
noise_sig = prop.sigma_sky

psf = galsim.Moffat(beta=3.5, scale_radius=0.8, flux=1.0, trunc=3)
psf_img = galsim.ImageD(stamp_size, stamp_size)
psf.drawImage(image=psf_img, scale=pixel_scale)
psf_ps = fq.pow_spec(psf_img.array)

if rank == 0:
    hdu = fits.PrimaryHDU(psf_img.array)
    psf_path = total_path + 'psf.fits'
    hdu.writeto(psf_path, overwrite=True)
    logger.info("size: %d, pixel_scale: %.2f, noise_sig: %.2f, total galaxy number: %d"%(stamp_size,
                                                                                         pixel_scale, noise_sig,
                                                                                         chips_num*divmod(cpus,shear_num)[0]))
logger.info("seed: %d"%seed)
fq = Fourier_Quad(stamp_size, 123)
t = 0
data_mea = numpy.zeros((chips_num*10000, 5))
for i in range(chips_num):
    t1 = time.clock()
    chip_path = total_path + "/chip/gal_chip_%s.fits"%(str(i+chip_s_id*chips_num).zfill(4))
    gal_pool = []
    logger.info("Start the %04d's chip..."%i)

    for k in range(10000):
        para_n = t+chip_s_id*chips_num*10000
        e1 = e1s[para_n]
        e2 = e2s[para_n]
        gal_flux = flux[para_n]
        ra = radius[para_n]
        btr = fbt[para_n][0]

        c_profile = numpy.random.randint(0, 10, 1)[0]
        if c_profile == 0:
            gal = galsim.DeVaucouleurs(half_light_radius=ra, trunc=4.5*ra).shear(e1=e1, e2=e2)
        else:
            bulge = galsim.Sersic(half_light_radius=0.6*ra, n=4, trunc=4.5*ra)# be careful
            disk = galsim.Sersic(half_light_radius=ra, n=1, trunc=4.5*ra)# be careful
            gal = bulge * btr + disk * (1-btr)
            gal = gal.shear(e1=e1, e2=e2)

        gal_s = gal.withFlux(gal_flux)
        gal_g = gal_s.shear(g1=g1, g2=g2)
        gal_c = galsim.Convolve([gal_g, psf])
        img = galsim.ImageD(stamp_size, stamp_size)
        gal_c.drawImage(image=img, scale=pixel_scale)
        gal_img = img.array + fq.draw_noise(0, noise_sig)
        gal_pool.append(gal_img)

        noise_img = fq.draw_noise(0,noise_sig)
        res = fq.shear_est(gal_img,psf_ps,noise_img,True)
        data_mea[t] = res[0], res[1], res[2], res[3], fq.flux2

        t += 1

    big_chip = fq.stack(gal_pool, 100)
    big_chip = numpy.float32(big_chip)
    hdu = fits.PrimaryHDU(big_chip)
    hdu.writeto(chip_path, overwrite=True)
    t2 = time.clock()
    logger.info("Finish the %04d's chip in %.2f sec"%(i, t2-t1))

if rank == 0:
    recv_buffer = numpy.empty((chips_num*10000*cpus, 5))
else:
    recv_buffer = None
comm.Gather(data_mea,recv_buffer, root=0)

if rank == 0:
    numpy.savez(result_path+"/data/data.npz",recv_buffer)
    FG1 = recv_buffer[:, 0]
    FG2 = recv_buffer[:, 1]
    FN = recv_buffer[:, 2]
    FU = recv_buffer[:, 3]
    DE1 = FN + FU
    DE2 = FN - FU
    g1_h, g1_h_sig = fq.fmin_g_new(FG1, DE1, bin_num=8)
    g2_h, g2_h_sig = fq.fmin_g_new(FG2, DE2, bin_num=8)
    print("g1: %.6f, mg1: %.6f, sig: %.6f, g2: %.6f, mg2: %.6f, sig: %.6f"%(g1, g1_h, g1_h_sig, g2, g2_h, g2_h_sig))
# ori_snr_path = result_path + "data/input_snr_%d.npz"%rank

te = time.clock()
logger.info("Used %.2f sec"%(te-ts))
