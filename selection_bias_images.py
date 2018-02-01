from sys import path
path.append('/home/hklee/work/fourier_quad/')
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import lsstetc
import time
from mpi4py import MPI
import h5py
import logging


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

with open("/home/hklee/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "total_data" in path:
        total_path = path.split("=")[1]
    elif "result" in path:
        result_path = path.split("=")[1]
    elif "parameter" in path:
        para_path = path.split("=")[1]
    elif "log" in path:
        log_path = path.split("=")[1]

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logfile = log_path + '%d_log.dat' %rank

lf = logging.FileHandler(logfile, 'w')
form = logging.Formatter('%(asctime)s - %(message)s')
lf.setFormatter(form)
logger.addHandler(lf)

stamp_size = 96
pixel_scale = 0.2
chips_num = 100
seed = rank*3424 + 53412
chip_s_id, shear_id = divmod(rank, 14)

fq = Fourier_Quad(stamp_size, seed)

shear_cata = para_path + "shear.dat"
shear = numpy.load(shear_cata)
g1 = shear["arr_0"][shear_id]
g2 = shear["arr_1"][shear_id]

paras = para_path + "para_%d.hdf5"%shear_id
f = h5py.File(paras,'r')
e1s = f["/e1"].value
e2s = f["/e2"].value
radius = f["/radius"].value
flux = f["/flux"].value
f.close()

prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
noise_sig = prop.sigma_sky

psf = galsim.Moffat(beta=3.5, scale_radius=1, flux=1.0)
psf_img = galsim.ImageD(stamp_size, stamp_size, pixel_scale)
psf.drawImage(image=psf_img)

#psf = psf_o.shear(e1=0.081, e2=-0.066)

if rank == 0:
    hdu = fits.PrimaryHDU(psf_img.array)
    psf_path = total_path + 'psf.fits'
    hdu.writeto(psf_path, overwrite=True)
    logger.info("size: %d, pixel_scale: %.2f, noise_sig: %.2f, total galaxy number: %d"%(stamp_size, pixel_scale, noise_sig, chips_num*divmod(cpus,14)[0]))
logger.info("seed: %d"%seed)

t = 0
for i in range(chips_num):
    t1 = time.time()
    chip_path = total_path + str(shear_id) + "/gal_chip_%s.fits"%(str(i+chip_s_id*chips_num).zfill(4))
    gal_pool = []
    logger.info("Start the %4d's chip..."%i)

    for k in range(10000):
        e1 = e1s[t+chip_s_id*chips_num]
        e2 = e2s[t+chip_s_id*chips_num]
        gal_flux = flux[t+chip_s_id*chips_num]
        ra = radius[t+chip_s_id*chips_num]
        t += 1

        bulge = galsim.Sersic(half_light_radius=ra-0.5, n=3.5)# be careful
        disk = galsim.Sersic(half_light_radius=ra, n=1.5)# be careful
        gal = bulge * 0.3 + disk * 0.7
        gal = gal.withFlux(gal_flux)

        gal_s = gal.shear(e1=e1, e2=e2)
        gal_g = gal_s.shear(g1=g1, g2=g2)
        gal_c = galsim.Convolve([psf, gal_g])
        img = galsim.ImageD(stamp_size, stamp_size, pixel_scale)
        gal_c.drawImage(image=img)
        gal_img = img.array + fq.draw_noise(0, noise_sig)
        gal_pool.append(gal_img)

    big_chip = fq.stack(gal_pool, 100)
    hdu = fits.PrimaryHDU(big_chip)
    hdu.writeto(chip_path, overwrite=True)
    t2 = time.clock()
    logger.info("Finish the %4d's chip in %.2f"%(i, t2-t1))

te = time.clock()
logger.info("Used %.2f sec"%(te-ts))
