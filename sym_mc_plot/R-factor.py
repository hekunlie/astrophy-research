import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
path.append('E:/Github/astrophy-research/my_lib/')
import time
from mpi4py import MPI
import tool_box
import numpy
from sys import argv
import galsim
from Fourier_Quad import Fourier_Quad
from astropy.io import fits

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

source = argv[1]

ini_path = "%s/work/envs/envs.dat"%my_home
total_path, para_path = tool_box.config(ini_path, ['get', 'get'],[['selection_bias', "%s_path"%source, '1'],
                                                    ['selection_bias', "%s_path_para"%source, '1']])

log_name = total_path + "logs/R_%02d.dat"%rank
logger = tool_box.get_logger(log_name)

# the parameters
para_contents = [["para","total_num",1], ["para","stamp_size",1], ["para", "stamp_col", 1], ["para","shear_num",1],
                 ["para","noise_sig",1], ["para", "pixel_scale", 1]]
para_items = tool_box.config(para_path+"para.ini", ['get', 'get', 'get', 'get', 'get', 'get'], para_contents)

total_chip_num = int(para_items[0])
stamp_size = int(para_items[1])
columns = int(para_items[2])
shear_num = int(para_items[3])
stamp_num = 10000

psf_path = total_path + "psf.fits"

chip_labels = tool_box.allot([i for i in range(total_chip_num)], cpus)[rank]
chip_num = len(chip_labels)

fq = Fourier_Quad(stamp_size, 123)
psf_img = galsim.Image(fits.open(psf_path)[0].data)
log_informs = "RANK: %d, SOURCE: %s, TOTAL CHIPS: %d, MY CHIPS: %d (%d ~ %d)"%(rank, source, total_chip_num,
                                                                               chip_num, chip_labels[0], chip_labels[-1])
logger.info(log_informs)
ts = time.time()
for i in range(shear_num):
    R_factor = numpy.zeros((chip_num*stamp_num, 1))
    for tag, j in enumerate(chip_labels):
        log_informs = "%02d/gal_chip_%04d.fits start"%(i, j)
        logger.info(log_informs)
        t1 = time.time()
        chip_path = total_path + "%d/gal_chip_%04d.fits"%(i, j)
        chip_img = fits.open(chip_path)[0].data
        gals = fq.segment(chip_img)
        for k in range(len(gals)):
            gal_img = galsim.Image(gals[k])
            R_factor[tag*stamp_num + k] = galsim.hsm.EstimateShear(gal_img, psf_img, strict=False).resolution_factor
        t2 = time.time()
        log_informs = "%02d/gal_chip_%04d.fits finish in %.2f"%(i, j, t2 - t1)
        logger.info(log_informs)
    if rank == 0:
        Recv_buffer = numpy.empty((total_chip_num*stamp_num, 1), dtype=numpy.float64)
    else:
        Recv_buffer = None
    comm.Gather(R_factor, Recv_buffer, root=0)
    if rank == 0:
        R_factor_path = total_path + "result/data/resolution_factor_%d.npz"%i
        numpy.savez(R_factor_path, Recv_buffer)

te = time.time()

log_informs = "RANK: %d, END: %.2f" % (rank, te - ts)
logger.info(log_informs)

