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

warnings.filterwarnings("ignore")
t1 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd, source, sex_filter = argv[1], argv[2], argv[3]

ini_path = "%s/work/envs/envs.dat"%my_home
log_name = "./m2_log.dat"
total_path, para_path = tool_box.config(ini_path, ['get', 'get'],
                                                   [['selection_bias', "%s_path"%source, '1'],
                                                    ['selection_bias', "%s_path_para"%source, '1']])

chip_num = 500
columns = 100
area_thresh = 6
max_distance = 8
size = int(tool_box.config(para_path+"para.ini", ["get"], [["para","size","1"]])[0])

if rank == 0:
    print(size)

total_fits_path = [total_path + "%d/gal_chip_%04d.fits"%(i, j) for i in range(14) for j in range(chip_num)]
total_cat_paths = [total_path + "result/data/%s/cat/%d_gal_chip_%04d.fits.cat"%(sex_filter,i, j) for i in range(14) for j in range(chip_num)]
allot_fits_path = tool_box.allot(total_fits_path, cpus)[rank]
allot_cat_path = tool_box.allot(total_cat_paths, cpus)[rank]

if cmd == "snr":
    for ii, chip_path in enumerate(allot_fits_path):
        cat_path = allot_cat_path[ii]
        # remove the previous cat file
        if os.path.exists(cat_path):
            os.remove(cat_path)
        comm = 'sex %s -CATALOG_NAME %s ' % (chip_path, cat_path)
        a = Popen(comm, shell=True)
        a.wait()

if cmd == "add" and rank < 14:
    snr_data_path = total_path + "result/data/%s/sex_%d.npz" %(sex_filter, rank)
    snr_data = numpy.zeros((chip_num * 10000, 6))
    for i in range(chip_num):
        cat_path = total_path + "result/data/%s/cat/%d_gal_chip_%04d.fits.cat"%(sex_filter,rank, i)
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

            if area >= area_thresh and distance <= max_distance:
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
    missing = ""
    for i in range(14):
        snr_data_path = total_path + "result/data/%s/sex_%d.npz" % (sex_filter, i)
        if not os.path.exists(snr_data_path):
            missing += " %d"%i
    log = "source: %s, filter: %s, time: %.2f, missing: %s"%(total_path.split("/")[-2], sex_filter, t2-t1, missing)
    tool_box.write_log(log_name, log)



