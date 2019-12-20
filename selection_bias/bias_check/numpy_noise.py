import numpy
from sys import path, argv
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
from mpi4py import MPI
import h5py
import tool_box
import time
from astropy.io import fits

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

total_path = argv[1]
log_path = total_path + "/logs/py_%02d.dat"%rank

stamp_size = 64
stamp_num = 10000
stamp_nx = 100
noise_sigma = 60
total_chips = 1000
shear_num = 10

logger = tool_box.get_logger(log_path)

m,n = divmod(total_chips, numprocs)
seed_step = 2
seed_gap = 2*seed_step*shear_num*(m+1)
seed = seed_gap*rank + 1

chips_label = numpy.zeros((numprocs, ),dtype=numpy.intc)+m
for i in range(n):
    chips_label[i] +=1

chip_st = chips_label[:rank].sum()
chip_ed = chip_st + chips_label[rank]

for i in range(shear_num):
    for j in range(chip_st, chip_ed):
        t1 = time.time()
        rng = numpy.random.RandomState(seed)

        noise1 = rng.normal(0, noise_sigma, (stamp_size*stamp_nx, stamp_size*stamp_nx))
        noise2 = rng.normal(0, noise_sigma, (stamp_size*stamp_nx, stamp_size*stamp_nx))

        # chip_path = total_path + "/%d/gal_chip_%04d.hdf5" % (i, j)
        # h5f = h5py.File(chip_path,"w")
        # h5f["/noise_1"] = noise1
        # h5f["/noise_2"] = noise2
        # h5f.close()

        chip_path = total_path + "/%d/gal_chip_%04d_0.fits" % (i, j)
        hdu = fits.PrimaryHDU(noise1)
        hdu.writeto(chip_path,overwrite=True)
        chip_path = total_path + "/%d/gal_chip_%04d_1.fits" % (i, j)
        hdu = fits.PrimaryHDU(noise2)
        hdu.writeto(chip_path,overwrite=True)

        t2 = time.time()
        logger.info("Shear %02d -- chip %04d, seed %d, %.2f sec"%(i, j, seed, t2-t1))
        seed += 1


