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
max_distance = 5
gal_num = 10000
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
    snr_data = numpy.zeros((chip_num * gal_num, 8))
    for i in range(chip_num):
        cat_path = total_path + "result/data/%s/cat/%d_gal_chip_%04d.fits.cat"%(sex_filter, rank, i)
        cata_data = numpy.loadtxt(cat_path)
        sex_data = tool_box.back_to_block(cata_data, gal_num, columns, size, size, 7, 6, 5, max_distance)
        snr_data[i * gal_num: (i + 1) * gal_num] = sex_data

    numpy.savez(snr_data_path, snr_data)
t2 = time.time()
if rank == 0:

    log = "source: %s, operation: %s filter: %s, time: %.2f, code: %s"%(total_path.split("/")[-2], argv[1], sex_filter, t2-t1, argv[0])
    tool_box.write_log(log_name, log)



