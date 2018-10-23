import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
path.append('E:/Github/astrophy-research/my_lib/')
import os
import time
from mpi4py import MPI
import tool_box
import numpy
import warnings
from sys import argv
import h5py
import matplotlib.pyplot as plt
import galsim
from Fourier_Quad import Fourier_Quad
from astropy.io import fits

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

source = argv[1]

total_chip_num = 500
columns = 100
gal_num = 10000
shear_num = 14

ini_path = "%s/work/envs/envs.dat"%my_home
log_name = "./r_log_%d.dat"%rank
total_path, para_path = tool_box.config(ini_path, ['get', 'get'],
                                                   [['selection_bias', "%s_path"%source, '1'],
                                                    ['selection_bias', "%s_path_para"%source, '1']])
psf_path = total_path + "psf.fits"
size = int(tool_box.config(para_path+"para.ini", ["get"], [["para","size","1"]])[0])

chip_labels = tool_box.allot(list(range(0, total_chip_num)), cpus)[rank]
chip_num = len(chip_labels)

fq = Fourier_Quad(size, 123)
psf_img = galsim.Image(fits.open(psf_path)[0].data)
log_informs = "RANK: %d, SOURCE: %s, TOTAL CHIPS: %d, MY CHIPS: %d (%d ~ %d)"%(rank, source, total_chip_num,
                                                                               chip_num, chip_labels[0], chip_labels[-1])

tool_box.write_log(log_name, log_informs)
ts = time.time()
for i in range(shear_num):
    R_factor = numpy.zeros((chip_num*gal_num, 1))
    for j in range(chip_num):
        log_informs = "%d / gal_chip_%04d.fits start"%(i, j)
        tool_box.write_log(log_name, log_informs)
        t1 = time.time()
        chip_path = total_path + "%d/gal_chip_%04d.fits"%(i, j)
        chip_img = fits.open(chip_path)[0].data
        gals = fq.segment(chip_img)
        for k in range(len(gals)):
            gal_img = galsim.Image(gals[k])
            result = galsim.hsm.EstimateShear(gal_img, psf_img, strict=False).resolution_factor
        t2 = time.time()
        log_informs = "%d / gal_chip_%04d.fits finish in %.2f"%(i, j, t2 - t1)
        tool_box.write_log(log_name, log_informs)
    if rank == 0:
        Recv_buffer = numpy.empty((total_chip_num*gal_num, 1), dtype=numpy.float64)
    else:
        Recv_buffer = None
    comm.Gather(R_factor, Recv_buffer, root=0)
    R_factor_path = total_path + "result/data/R_fractor_%d.npz"%i
    numpy.savez(R_factor_path, Recv_buffer)

te = time.time()

log_informs = "RANK: %d, END: %.2f" % (rank, te - ts)
tool_box.write_log(log_name, log_informs)

