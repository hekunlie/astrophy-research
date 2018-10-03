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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ini_path = "%s/work/envs/envs.dat"%my_home
total_path = tool_box.config(ini_path,['get'], [['selection_bias', "pts_path", '1']])[0]
fits_path = total_path + "%d/"%rank

chip_num = 500
data = numpy.zeros((chip_num*10000, 6))
columns = 100
for i in range(chip_num):
    # SNR estimation
    chip_path = fits_path + "gal_chip_%04d"%i
    cat_path = chip_path + ".cat"

    # remove the previous cat file
    if os.path.exists(cat_path):
        os.remove(cat_path)

    comm = 'sex %s -CATALOG_NAME %s.cat '%(chip_path, cat_path)
    a = Popen(comm, shell=True)
    a.wait()

    # add to catalog
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
    for i in range(len(cata_data)):
        xx = a_x[i] - 1
        yy = a_y[i] - 1
        snr = a_snr[i]
        area = a_area[i]

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
    data[k * 10000:(k + 1) * 10000] = sex_data




