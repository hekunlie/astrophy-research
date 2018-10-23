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
import h5py
import matplotlib.pyplot as plt

warnings.filterwarnings("ignore")
t1 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd, source, sex_filter, max_distance = argv[1], argv[2], argv[3], float(argv[4])

ini_path = "%s/work/envs/envs.dat"%my_home
log_name = "./m2_log.dat"
total_path, para_path = tool_box.config(ini_path, ['get', 'get'],
                                                   [['selection_bias', "%s_path"%source, '1'],
                                                    ['selection_bias', "%s_path_para"%source, '1']])

chip_num = 500
columns = 100
area_thresh = 6
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

if cmd == "check" and rank < 14:
    # check
    check_path = total_path + "data/result/check/%s/"%sex_filter
    if os.path.exists(check_path):
        os.rmdir(check_path)
    else:
        os.makedirs(check_path)

    input_para_path = para_path + "para_%d.hdf5"%rank
    para_h5 = h5py.File(input_para_path, "r")
    input_mag = para_h5["/mag"].value
    input_flux = para_h5["/flux"].value
    para_h5.close()

    snr_data_path = total_path + "result/data/%s/sex_%d.npz" %(sex_filter, rank)
    snr_data = numpy.load(snr_data_path)["arr_0"]
    sex_snr = snr_data[:, 0]
    sex_mag = snr_data[:, 2]
    sex_area = snr_data[:, 5]
    idx = sex_mag > 0

    plt.figure(figsize=(21, 21))
    plt.subplot(331)
    plt.scatter(sex_snr[idx], input_mag[idx], s=1)
    plt.xlabel("SNR")
    plt.ylabel("INPUT MAG")
    plt.subplot(332)
    plt.scatter(sex_mag[idx], input_mag[idx], s=1)
    plt.xlabel("MAG_AUTO")
    plt.ylabel("INPUT MAG")
    plt.subplot(333)
    plt.scatter(sex_area[idx], input_mag[idx], s=1)
    plt.xlabel("SEX_AREA")
    plt.ylabel("INPUT MAG")
    plt.subplot(334)
    idx_s = sex_snr < 200
    plt.hist(sex_snr[idx&idx_s], 100)
    plt.xlabel("SEX_SNR")
    plt.subplot(335)
    plt.hist(sex_mag[idx], 100)
    plt.xlabel("MAG_AUTO")
    plt.subplot(336)
    plt.hist(input_mag, 100)
    plt.xlabel("INPUT MAG")

    sig_level = sex_filter.split("_")[-1]
    meas_para_path = total_path + "result/data/data_%ssig/data_%d_0.hdf5" % (sig_level, rank)
    if os.path.exists(meas_para_path):
        meas_para = h5py.File(meas_para_path, "r")
        fsnr = meas_para["/data"].value[:, 3]
        print(fsnr.shape, input_mag.shape)
        meas_para.close()
        plt.subplot(337)
        plt.scatter(fsnr, input_mag, s=1)
        plt.xlabel("F-SNR (1.5$/sigma$)")
        plt.ylabel("INPUT MAG")
        plt.hist(fsnr, 100)
        plt.xlabel("F-SNR (%s$/sigma$)"%sig_level)
        plt.subplot(339)
        plt.scatter(sex_snr, fsnr, s=1)
        plt.xlabel("F-SNR (%s$/sigma$)"%sig_level)
        plt.ylabel("SEX_SNR (%s$/sigma$)"%sig_level)
    pic_nm = check_path + "check_%d.png"%rank
    plt.savefig(pic_nm)
    plt.close()
t2 = time.time()
if rank == 0:
    log = "operation: %s, source: %s, max radius: %.2f, stamp size: %d, filter: %s, time: %.2f, code: %s"\
          %(argv[1], total_path.split("/")[-2], max_distance, size, sex_filter, t2-t1, argv[0])
    tool_box.write_log(log_name, log)



