import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/fourier_quad/'%my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import time
from mpi4py import MPI
import h5py
import tool_box
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()
source = argv[1]
envs_path = "%s/work/envs/envs.dat"%my_home
get_contents = [['selection_bias', "%s_path"%source, '1'],['selection_bias', "%s_path_result"%source, '1'],
                ['selection_bias', "%s_path_para"%source, '1'],['selection_bias', "%s_path_log"%source, '1']]
path_items = tool_box.config(envs_path,['get','get','get','get'], get_contents)
total_path, result_path, para_path, log_path = path_items

logger = tool_box.get_logger(log_path+"%d_logs.dat"%rank)

stamp_size = 100
stamp_col = 100
pixel_scale = 0.187
shear_num = 10
noise_sig = 65
total_chips_num = 20
seed = rank*344 + 1212

ny, nx = stamp_col*stamp_size, stamp_col*stamp_size
fq = Fourier_Quad(stamp_size, seed)

psf = galsim.Moffat(beta=3.5, fwhm=1, flux=1.0, trunc=2)
#psf = psf_o.shear(e1=0.081, e2=-0.066)

if rank == 0:
    psf_img = galsim.ImageD(stamp_size, stamp_size)
    psf.drawImage(image=psf_img, scale=pixel_scale)
    hdu = fits.PrimaryHDU(psf_img.array)
    psf_path = total_path + 'psf.fits'
    hdu.writeto(psf_path, overwrite=True)
    logger.info("desti: %s, size: %d, pixel_scale: %.3f, noise_sig: %.2f, total galaxy number: %d"
                %(source,stamp_size, pixel_scale, noise_sig, total_chips_num))
logger.info("seed: %d"%seed)

chip_tags = [i for i in range(total_chips_num)]
chip_tags_rank = tool_box.allot(chip_tags, cpus)[rank]

for shear_id in range(shear_num):
    shear_cata = para_path + "shear.npz"
    shear = numpy.load(shear_cata)
    g1 = shear["arr_0"][shear_id]
    g2 = shear["arr_1"][shear_id]
    if rank == 0:
        print(g1,g2, type(g1))
    paras = para_path + "para_%d.hdf5" % shear_id
    f = h5py.File(paras, 'r')
    e1s = f["/e1"].value
    e2s = f["/e2"].value
    radius = f["/radius"].value
    flux = f["/flux"].value
    fbt = f['/btr'].value
    gal_profile = f["/type"].value
    f.close()
    #
    rng = numpy.random.RandomState(seed + shear_id)
    # t = 1
    # fig = plt.figure(figsize=(int(5*len(chip_tags_rank)), 5))
    for chip_tag in chip_tags_rank:
        t1 = time.clock()
        chip_path = total_path + "%s/gal_chip_%04d.fits"%(shear_id, chip_tag)
        gal_pool = []
        para_n = chip_tag*10000
        logger.info("Start the %04d's chip at shear %02d...\n"
                    "Para from: %d" % (chip_tag, shear_id, para_n))
        # ax = fig.add_subplot(1, int(len(chip_tags_rank)), t)
        # t+=1
        para_no = []
        prof_no = []
        for k in range(10000):
            e1 = e1s[para_n+k, 0]
            e2 = e2s[para_n+k, 0]
            gal_flux = flux[para_n+k, 0]
            ra = radius[para_n+k, 0]
            btr = fbt[para_n+k, 0]
            c_profile = gal_profile[para_n+k, 0]
            para_no.append(para_n+k)
            prof_no.append(c_profile)

            if c_profile == 1:
                gal = galsim.DeVaucouleurs(half_light_radius=ra, trunc=4.5*ra, flux_untruncated=False).shear(e1=e1, e2=e2)
            else:
                bulge = galsim.Sersic(half_light_radius=ra, n=4, trunc=4.5*ra, flux_untruncated=False)# be careful
                disk = galsim.Sersic(scale_radius=ra, n=1, trunc=4.5*ra, flux_untruncated=False)# be careful
                gal = bulge * btr + disk * (1-btr)
                gal = gal.shear(e1=e1, e2=e2)

            gal_s = gal.withFlux(gal_flux)
            gal_g = gal_s.shear(g1=g1, g2=g2)
            gal_c = galsim.Convolve([gal_g, psf])

            img = galsim.ImageD(stamp_size, stamp_size)
            gal_c.drawImage(image=img, scale=pixel_scale)
            gal_pool.append(img.array)

        noise_img = rng.normal(0, noise_sig, nx*ny).reshape((ny, nx))
        big_chip = fq.stack(gal_pool, stamp_col) + noise_img
        big_chip = numpy.float32(big_chip)
        hdu = fits.PrimaryHDU(big_chip)
        hdu.writeto(chip_path, overwrite=True)
        t2 = time.clock()
        logger.info("Finish the %04d's chip at shear %02d in %.2f sec"%(chip_tag, shear_id, t2-t1))
        # ax.scatter(para_no, prof_no)
    # pic_path = "/home/hkli/work/test/pic/%d_%d.png"%(shear_id, rank)
    # plt.savefig(pic_path)
    # plt.close()
te = time.clock()
logger.info("Used %.2f sec. %s"%(te-ts, total_path))
