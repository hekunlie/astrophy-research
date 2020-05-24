import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import time
from mpi4py import MPI
import h5py
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

parent_path = argv[1]
tag = int(argv[2])

gal_fluxs = [3000, 6000, 9000, 16000, 800000]

seed_step = 1

sersic_idx = 0.31
scale_radius = 0.3
gal_flux = gal_fluxs[tag]
noise_sig = 60
shear_num = 40
total_chips = 4000
pixel_scale = 0.187
stamp_size = 48
stamp_nx = 100
stamp_ny = 100
stamp_num = stamp_nx * stamp_ny

nx, ny = stamp_nx*stamp_size, stamp_ny*stamp_size

fq = Fourier_Quad(stamp_size, 123)

logger = tool_box.get_logger(parent_path+"/logs/%d_logs.dat" % rank)

psf = galsim.Moffat(beta=3.5, fwhm=0.7, flux=1.0, trunc=1.4).shear(e1=0.1, e2=0.)
if rank == 0:
    psf_img = galsim.ImageD(stamp_size, stamp_size)
    psf.drawImage(image=psf_img, scale=pixel_scale)
    psf_arr = numpy.float32(psf_img.array)
    hdu = fits.PrimaryHDU(psf_arr)
    psf_path = parent_path + '/psf.fits'
    hdu.writeto(psf_path, overwrite=True)


gal = galsim.Sersic(scale_radius=scale_radius, n=sersic_idx, trunc=4.5 * scale_radius, flux=1.0)

# task allocation
chip_tags = [i for i in range(total_chips)]
chip_tags_rank = tool_box.alloc(chip_tags, cpus)[rank]

shear_path = parent_path + "/parameters/shear.hdf5"
h5f = h5py.File(shear_path, "r")
g1_input = h5f["/g1"][()]
g2_input = h5f["/g2"][()]
h5f.close()

img_buffer = numpy.zeros((ny, nx))
counts = 0

for shear_id in range(shear_num):

    if rank == 0:
        if not os.path.exists(parent_path + "/imgs/%d"%shear_id):
            os.makedirs(parent_path + "/imgs/%d"%shear_id)
    comm.Barrier()

    g1 = g1_input[shear_id]
    g2 = g2_input[shear_id]

    paras = parent_path + "/parameters/para_%d.hdf5" % shear_id
    h5f = h5py.File(paras, 'r')
    e1s = h5f["/e1"][()]
    e2s = h5f["/e2"][()]
    h5f.close()

    logger.info("SHEAR ID: %02d, RANk: %02d, e1s: %.3f, e2s: %.3f"
                %(shear_id, rank, e1s[int(0.5*total_chips)], e2s[int(0.9*total_chips)]))

    for t, chip_tag in enumerate(chip_tags_rank):
        t1 = time.clock()

        chip_path = parent_path + "/imgs/%d/chip_%04d.fits" % (shear_id, chip_tag)

        logger.info("SHEAR ID: %02d, Start the %04d's chip." % (shear_id, chip_tag))

        para_n = chip_tag * stamp_num

        for k in range(stamp_num):
            e1 = e1s[para_n + k]
            e2 = e2s[para_n + k]

            gal_e = gal.shear(e1=e1, e2=e2).withFlux(gal_flux)
            gal_s = gal_e.shear(g1=g1, g2=g2)
            gal_c = galsim.Convolve([gal_s, psf])
            gal_img = galsim.ImageD(stamp_size, stamp_size)
            gal_c.drawImage(image=gal_img, scale=pixel_scale)
            iy, ix = divmod(k, stamp_nx)
            fq.stack_new(img_buffer, gal_img.array, iy, ix)

        big_chip = numpy.float32(img_buffer)
        hdu = fits.PrimaryHDU(big_chip)
        hdu.writeto(chip_path, overwrite=True)
        t2 = time.clock()
        logger.info("SHEAR ID: %02d, Finish the %04d's chip in %.2f sec" % (shear_id, chip_tag, t2 - t1))

logger.info("Flux: %.1f"%gal_flux)

