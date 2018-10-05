import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
from subprocess import Popen
import os
import time
from mpi4py import MPI
import tool_box
import numpy
import warnings
from sys import argv
import datetime


warnings.filterwarnings("ignore")
t1 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]

ini_path = "%s/work/envs/envs.dat"%my_home
source = "dimmerm2"
total_path, para_path = tool_box.config(ini_path, ['get', 'get'],
                                                   [['selection_bias', "%s_path"%source, '1'],
                                                    ['selection_bias', "%s_path_para"%source, '1']])

chip_num = 500
columns = 100
size = int(tool_box.config(para_path+"para.ini", ["get"], [["para","size","1"]])[0])

if rank == 0:
    print(size)

fits_paths = [total_path + "%d/gal_chip_%04d.fits"%(i, j) for i in range(14) for j in range(chip_num)]
fits_path = tool_box.allot(fits_paths, cpus)[rank]

if cmd == "snr":
    for chip_path in fits_path:
        cat_path = chip_path + ".cat"
        # remove the previous cat file
        if os.path.exists(cat_path):
            os.remove(cat_path)
        comm = 'sex %s -CATALOG_NAME %s ' % (chip_path, cat_path)
        a = Popen(comm, shell=True)
        a.wait()

if rank < 14:
    snr_data_path = total_path + "result/data/sex3_2_%d.npz" %rank
    snr_data = numpy.zeros((chip_num * 10000, 6))
    for i in range(chip_num):
        cat_path = total_path + "%d/gal_chip_%04d.fits.cat"%(rank, i)
        cata_data = numpy.loadtxt(cat_path)
        a_snr = cata_data[:, 0]
        a_mag_iso = cata_data[:, 1]
        a_mag_auto = cata_data[:, 2]
        a_mag_petro = cata_data[:, 3]
        a_mag_win = cata_data[:, 4]
        a_area = cata_data[:, 5]
        a_x = cata_data[:, 6]
        a_y = cata_data[:, 7]
        sex_data = numpy.zeros((10000, 6))
        for ii in range(len(cata_data)):
            xx = a_x[ii] - 1
            yy = a_y[ii] - 1
            snr = a_snr[ii]
            area = a_area[ii]

            mx, modx = divmod(xx, size)
            my, mody = divmod(yy, size)
            tag = int(columns * my + mx)
            distance = numpy.sqrt((modx - size / 2) ** 2 + (mody - size / 2) ** 2)

            if area >= 6 and distance <= 8:
                if area > sex_data[tag, 1]:
                    sex_data[tag, 0] = snr
                    sex_data[tag, 1] = area
                    sex_data[tag, 2] = a_mag_iso[i]
                    sex_data[tag, 3] = a_mag_auto[i]
                    sex_data[tag, 4] = a_mag_petro[i]
                    sex_data[tag, 5] = a_mag_win[i]
        snr_data[i * 10000:(i + 1) * 10000] = sex_data

    numpy.savez(snr_data_path, snr_data)
t2 = time.time()
if rank == 0:
    log = "%s, %s, %.2f"%( total_path.split("/")[-2], snr_data_path.split("/")[-1], t2-t1)
    tool_box.write_log("./m2_log.dat",log)



