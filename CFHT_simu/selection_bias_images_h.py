import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import time
from mpi4py import MPI
import h5py
import tool_box
import shutil

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

source = argv[1]

envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['selection_bias', "%s_path" % source, '1'], ['selection_bias', "%s_path_result" % source, '1'],
                ['selection_bias', "%s_path_para" % source, '1'], ['selection_bias', "%s_path_log" % source, '1']]
path_items = tool_box.config(envs_path, ['get', 'get', 'get', 'get'], get_contents)
total_path, result_path, para_path, log_path = path_items

logger = tool_box.get_logger(log_path + "%d_logs.dat" % rank)

# the parameters
para_contents = [["para","total_num",1], ["para","stamp_size",1], ["para", "stamp_col", 1], ["para","shear_num",1],
                 ["para","noise_sig",1], ["para", "pixel_scale", 1]]
para_items = tool_box.config(para_path+"para.ini", ['get', 'get', 'get', 'get', 'get', 'get'], para_contents)

total_chips_num = int(para_items[0])
stamp_size = int(para_items[1])
stamp_col = int(para_items[2])
shear_num = int(para_items[3])
noise_sig = int(para_items[4])
pixel_scale = float(para_items[5])
stamp_num = 10000

finish_path = "%s/work/test/job/%s/finish_%d.dat"%(my_home, source, rank)
if rank == 0:
    indicator = "%s/work/test/job/%s"%(my_home, source)
    if os.path.exists(indicator):
        shutil.rmtree(indicator)
    os.makedirs(indicator)
comm.Barrier()

total_gal_num = total_chips_num * stamp_num
seed = rank * 344 + 11121

ny, nx = stamp_col * stamp_size, stamp_col * stamp_size
fq = Fourier_Quad(stamp_size, seed)

# PSF
psf = galsim.Moffat(beta=3.5, fwhm=0.7, flux=1.0, trunc=2)
if rank == 0:
    psf_img = galsim.ImageD(stamp_size, stamp_size)
    psf.drawImage(image=psf_img, scale=pixel_scale)
    hdu = fits.PrimaryHDU(psf_img.array)
    psf_path = total_path + 'psf.fits'
    hdu.writeto(psf_path, overwrite=True)
    logger.info("desti: %s, size: %d, pixel_scale: %.3f, noise_sig: %.2f, total galaxy number: %d"
                %(source,stamp_size, pixel_scale, noise_sig, total_chips_num))
logger.info("seed: %d"%seed)

# task allocation
chip_tags = [i for i in range(total_chips_num)]
chip_tags_rank = tool_box.allot(chip_tags, cpus)[rank]

for shear_id in range(shear_num):
    shear_cata = para_path + "shear.npz"
    shear = numpy.load(shear_cata)
    g1 = shear["arr_0"][shear_id]
    g2 = shear["arr_1"][shear_id]

    paras = para_path + "para_%d.hdf5" % shear_id
    f = h5py.File(paras, 'r')
    e1s = f["/e1"].value
    e2s = f["/e2"].value
    radius = f["/radius"].value
    flux = f["/flux"].value
    fbt = f['/btr'].value
    gal_profile = f["/type"].value

    # for checking
    logger.info("SHEAR ID: %02d, RANk: %02d, e1s: %.3f, e2s: %.3f, radius: %.2f, fbt: %.2f"
                %(shear_id, rank, e1s[int(0.5*total_gal_num)], e2s[int(0.9*total_gal_num)],
                  radius[int(0.1*total_gal_num)], fbt[int(0.99*total_gal_num)]))

    for t, chip_tag in enumerate(chip_tags_rank):
        t1 = time.clock()

        rng = numpy.random.RandomState(seed + shear_id + t)

        chip_path = total_path + "%s/gal_chip_%04d.fits" % (shear_id, chip_tag)
        gal_pool = []
        logger.info("SHEAR ID: %02d, Start the %04d's chip.." % (shear_id, chip_tag))

        para_n = chip_tag * stamp_num

        for k in range(stamp_num):
            e1 = e1s[para_n + k, 0]
            e2 = e2s[para_n + k, 0]
            gal_flux = flux[para_n + k, 0]
            ra = radius[para_n + k, 0]
            btr = fbt[para_n + k, 0]

            # regular galaxy
            c_profile = gal_profile[para_n + k, 0]
            if c_profile == 1:
                gal = galsim.DeVaucouleurs(half_light_radius=ra, trunc=4.5 * ra,flux=1.0)
            else:
                bulge = galsim.Sersic(half_light_radius=ra, n=4, trunc=4.5 * ra,flux=1.0)  # be careful
                disk = galsim.Sersic(scale_radius=ra, n=1, trunc=4.5 * ra,flux=1.0)  # be careful
                gal = bulge * btr + disk * (1 - btr)
            gal_e = gal.shear(e1=e1, e2=e2)
            gal_f = gal_e.withFlux(gal_flux)

            # random walk
            # gal_rng_seed = seed + shear_id + t + k
            # gal_rng = galsim.BaseDeviate(gal_rng_seed)
            # gal_f = galsim.randwalk.RandomWalk(npoints=80, half_light_radius=ra, flux=gal_flux, rng=gal_rng)#.shear(e1=e1, e2=e2)

            gal_s = gal_f.shear(g1=g1, g2=g2)
            gal_c = galsim.Convolve([gal_s, psf])

            img = galsim.ImageD(stamp_size, stamp_size)
            gal_c.drawImage(image=img, scale=pixel_scale)
            gal_pool.append(img.array)

        noise_img = rng.normal(0, noise_sig, nx * ny).reshape((ny, nx))
        big_chip = fq.stack(gal_pool, stamp_col) + noise_img
        big_chip = numpy.float32(big_chip)
        hdu = fits.PrimaryHDU(big_chip)
        hdu.writeto(chip_path, overwrite=True)
        t2 = time.clock()
        logger.info("SHEAR ID: %02d, Finish the %04d's chip in %.2f sec" % (shear_id, chip_tag, t2 - t1))


with open(finish_path, "w") as f:
    f.write("done")

te = time.clock()
logger.info("Used %.2f sec. %s" % (te - ts, total_path))
