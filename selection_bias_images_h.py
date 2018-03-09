from sys import path
path.append('/home/hkli/work/fourier_quad')
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

with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "total" in path:
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

stamp_size = 90
pixel_scale = 0.2
chips_num = int(500/int(cpus/14))
seed = rank*34424 + 41112
chip_s_id, shear_id = divmod(rank, 14)

fq = Fourier_Quad(stamp_size, seed)

shear_cata = para_path + "shear.npz"
shear = numpy.load(shear_cata)
g1 = shear["arr_0"][shear_id]
g2 = shear["arr_1"][shear_id]

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

#psf = psf_o.shear(e1=0.081, e2=-0.066)

if rank == 0:
    hdu = fits.PrimaryHDU(psf_img.array)
    psf_path = total_path + 'psf.fits'
    hdu.writeto(psf_path, overwrite=True)
    logger.info("size: %d, pixel_scale: %.2f, noise_sig: %.2f, total galaxy number: %d"%(stamp_size, pixel_scale, noise_sig, chips_num*divmod(cpus,14)[0]))
logger.info("seed: %d"%seed)

t = 0

for i in range(chips_num):
    t1 = time.clock()
    chip_path = total_path + str(shear_id) + "/gal_chip_%s.fits"%(str(i+chip_s_id*chips_num).zfill(4))
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
            gal = galsim.DeVaucouleurs(half_light_radius=ra, trunc=5*ra).shear(e1=e1, e2=e2)
        else:
            bulge = galsim.Sersic(half_light_radius=0.6*ra, n=4, trunc=5*ra)# be careful
            disk = galsim.Sersic(half_light_radius=ra, n=1, trunc=5*ra)# be careful
            gal = bulge * btr + disk * (1-btr)
            gal = gal.shear(e1=e1, e2=e2)
        t += 1

        gal_s = gal.withFlux(gal_flux)
        gal_g = gal_s.shear(g1=g1, g2=g2)
        gal_c = galsim.Convolve([gal_g, psf])
        img = galsim.ImageD(stamp_size, stamp_size)
        gal_c.drawImage(image=img, scale=pixel_scale)
        gal_img = img.array + fq.draw_noise(0, noise_sig)
        gal_pool.append(gal_img)

    big_chip = fq.stack(gal_pool, 100)
    # big_chip = numpy.float32(big_chip)
    hdu = fits.PrimaryHDU(big_chip)
    hdu.writeto(chip_path, overwrite=True)
    t2 = time.clock()
    logger.info("Finish the %04d's chip in %.2f sec"%(i, t2-t1))

te = time.clock()
logger.info("Used %.2f sec"%(te-ts))
