import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
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

cmd, total_path, sex_filter, max_distance = argv[1], argv[2], argv[3], float(argv[4])

log_name = "./m2_log.dat"
logger = tool_box.get_logger(log_name)
para_path = total_path + "/parameters"
# print(total_path, para_path)
para_contents = [["para", "total_num", 1], ["para", "stamp_size", 1],
                 ["para", "stamp_col", 1], ["para", "shear_num", 1],
                 ["para", "noise_sig", 1]]
para_items = tool_box.config(para_path+"/para.ini", ['get', 'get', 'get', 'get','get'], para_contents)

chip_num = int(para_items[0])
size = int(para_items[1])
columns = int(para_items[2])
shear_num = int(para_items[3])
noise_sig = float(para_items[4])

area_thresh = 5
gal_num = 10000

snr_idx = 0
flux_auto_idx = 1
flux_err_idx = 2
mag_auto_idx = 3
area_idx = 4
x_idx = 5
y_idx = 6

if rank == 0:
    log = "START: operation: %s, source: %s, code: %s, shear num: %d"\
          %(argv[1], total_path.split("/")[-1], argv[0], shear_num)
    logger.info(log)

if cmd == "snr":
    ##################################### measure SNR ######################################################
    # this section will call sextractor to detect the source in fits files

    filter_names = ["sex2_4", "sex4_4",
                    "sex2_2", "sex4_2",
                    "sex2_1.5", "sex4_1.5"]
    gauss_filters = ["gauss_2.0_5x5", "gauss_4.0_7x7",
                     "gauss_2.0_5x5", "gauss_4.0_7x7",
                     "gauss_2.0_5x5", "gauss_4.0_7x7"]

    if sex_filter not in filter_names:
        print("%s not found in "%sex_filter, filter_names)
        raise KeyError

    total_fits_path = [total_path + "/%d/gal_chip_%04d.fits" % (i, j) for i in range(shear_num) for j in
                       range(chip_num)]
    total_cat_paths = [total_path + "/result/data/%s/cat/%d_gal_chip_%04d.fits.cat" % (sex_filter, i, j)
                       for i in range(shear_num) for j in range(chip_num)]
    allot_fits_path = tool_box.allot(total_fits_path, cpus)[rank]
    allot_cat_path = tool_box.allot(total_cat_paths, cpus)[rank]

    if rank == 0:
        with open("./default.sex_ori", "r") as f:
            contents = f.readlines()
        sig_level = float(sex_filter.split("_")[1])*noise_sig
        contents[16] = "DETECT_THRESH    %.2f          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n"%sig_level
        contents[18] = "ANALYSIS_THRESH  %.2f          # <sigmas> or <threshold>,<ZP> in mag.arcsec-2\n"%sig_level
        contents[23] = "FILTER_NAME      %s.conv   # name of the file containing the filter\n"%gauss_filters[filter_names.index(sex_filter)]
        with open("./default.sex", "w") as f:
            f.writelines(contents)

        cat_path = total_path + "/result/data/%s/cat/" % sex_filter
        if not os.path.exists(cat_path):
            os.mkdir(cat_path)
        logger.info("%s, %s"%(sex_filter,gauss_filters[filter_names.index(sex_filter)]))
        # print(sex_filter,gauss_filters[filter_names.index(sex_filter)])
        # time.sleep(30)
    comm.Barrier()

    for ii, chip_path in enumerate(allot_fits_path):
        cat_path = allot_cat_path[ii]
        # remove the previous cat file
        if os.path.exists(cat_path):
            os.remove(cat_path)
        comm = "sex %s -CATALOG_NAME %s" % (chip_path, cat_path)
        a = Popen(comm, shell=True)
        a.wait()


##################################### add to catalog ######################################################
if cmd == "add" and rank < shear_num:
    # each rank only operate shear_id == rank !!!
    snr_data_path = total_path + "/result/data/%s/sex_%d.npz" %(sex_filter, rank)
    snr_data = numpy.zeros((chip_num * gal_num, 7))
    cent = size/2. - 0.5
    logger.info("%02d start..."%rank)
    for i in range(chip_num):
        cat_path = total_path + "/result/data/%s/cat/%d_gal_chip_%04d.fits.cat"%(sex_filter, rank, i)
        cata_data = numpy.loadtxt(cat_path)
        sex_data = tool_box.back_to_block(cata_data, gal_num, columns, size, size, cent, cent, y_idx, x_idx, area_idx, max_distance)
        snr_data[i * gal_num: (i + 1) * gal_num] = sex_data
    numpy.savez(snr_data_path, snr_data)
    logger.info("%02d Begin to write to disk..."%rank)

    # if the Pk0 data don't exist, measure them first
    fourier_path = "%s/result/data/data_1.5sig/data_%d.hdf5" % (total_path, rank)
    h5f = h5py.File(fourier_path, "r")
    fourier_data = h5f["/data"][()]
    h5f.close()

    # Pk0
    h5path = total_path + "/result/data/%s/flux2_ex1_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = fourier_data[:, 4]
    h5f.close()

    # Pk0_fit
    h5path = total_path + "/result/data/%s/flux2_ex2_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = fourier_data[:, 5]
    h5f.close()

    # max(Pk0,Pk0_fit)
    h5path = total_path + "/result/data/%s/flux2_ex3_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = fourier_data[:, 6]
    h5f.close()

    # mask
    mask = numpy.ones((chip_num * gal_num,), dtype=numpy.intc)
    idx = snr_data[:, snr_idx] <= 0
    mask[idx] = 0
    h5path = total_path + "/result/data/%s/mask_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = mask
    h5f.close()

    # true magnitude
    h5path = total_path + "/parameters/para_%d.hdf5" % rank
    h5f = h5py.File(h5path, "r")
    mag_true = h5f["/mag"][()]
    mag_true.shape = (mag_true.shape[0],)
    h5f.close()
    h5path = total_path + "/result/data/%s/mag_true_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = -mag_true
    h5f.close()

    # magnitude
    h5path = total_path + "/result/data/%s/mag_auto_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = -snr_data[:, mag_auto_idx]
    h5f.close()

    # snr
    h5path = total_path + "/result/data/%s/snr_sex_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = snr_data[:, snr_idx]
    h5f.close()

    # snr_auto
    snr_data[:, flux_err_idx][idx] = 1
    h5path = total_path + "/result/data/%s/snr_auto_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = snr_data[:, flux_auto_idx] / snr_data[:, flux_err_idx]
    h5f.close()

    # area
    h5path = total_path + "/result/data/%s/area_%d.hdf5" % (sex_filter, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = snr_data[:, area_idx]
    h5f.close()

    logger.info("%02d Begin to write to disk..." % rank)

##################################### check ######################################################
if cmd == "check" and rank < shear_num:
    # check
    check_path = total_path + "/result/data/check/%s/"%sex_filter
    if rank == 0:
        if not os.path.exists(check_path):
            # os.removedirs(check_path)
            os.makedirs(check_path)

    # it has been found that the program will be hung on until the bcast() has been done
    # even if the 'cmd = "done"' is putted before the 'if ...'
    # and all the threads will synchronized before the communion, the thread that runs faster will be hung on
    # so a judge of the 'cmd' is not need
    comm.Barrier()

    input_para_path = para_path + "/para_%d.hdf5"%rank
    para_h5 = h5py.File(input_para_path, "r")
    input_mag = para_h5["/mag"][()]
    input_flux = para_h5["/flux"][()]
    para_h5.close()

    snr_data_path = total_path + "/result/data/%s/sex_%d.npz" %(sex_filter, rank)
    snr_data = numpy.load(snr_data_path)["arr_0"]
    sex_snr = snr_data[:, snr_idx]
    sex_mag = snr_data[:, mag_auto_idx]
    sex_area = snr_data[:, area_idx]
    idx = sex_mag > 0
    ms = 0.2
    plt.figure(figsize=(35, 14))
    plt.subplot(2,5,1)
    plt.scatter(sex_snr[idx], input_mag[idx], s=ms)
    plt.xlabel("SNR")
    plt.ylabel("INPUT MAG")
    plt.title("%d detected"%len(sex_snr[idx]))
    plt.subplot(2,5,2)
    plt.scatter(sex_mag[idx], input_mag[idx], s=ms)
    plt.xlabel("MAG_AUTO")
    plt.ylabel("INPUT MAG")
    plt.subplot(2,5,3)
    plt.scatter(sex_area[idx], input_mag[idx], s=ms)
    plt.xlabel("SEX_AREA")
    plt.ylabel("INPUT MAG")
    plt.subplot(2,5,6)
    idx_s = sex_snr < 200
    plt.hist(sex_snr[idx&idx_s], 100)
    plt.xlabel("SEX_SNR")
    plt.subplot(2,5,7)
    plt.hist(sex_mag[idx], 100)
    plt.xlabel("MAG_AUTO")
    plt.subplot(2,5,8)
    plt.hist(input_mag, 100)
    plt.xlabel("INPUT MAG")

    sig_level = float(sex_filter.split("_")[1])
    meas_para_path = total_path + "/result/data/data_%.1fsig/data_%d.hdf5" % (sig_level, rank)
    if os.path.exists(meas_para_path):
        meas_para = h5py.File(meas_para_path, "r")
        f_data = meas_para["data"].value
        fsnr = f_data[:,3]
        detect = f_data[:,0]
        idx_f = detect > 0
        meas_para.close()

        plt.subplot(2,5,4)
        plt.scatter(fsnr[idx_f], input_mag[idx_f], s=ms)
        plt.xlabel("F-SNR (%s$\sigma$)"%sig_level)
        plt.ylabel("INPUT MAG")

        plt.subplot(2,5,9)
        idx_f_lim = fsnr < 30
        plt.hist(fsnr[idx_f_lim&idx_f], 100)
        plt.xlabel("F-SNR (%s$\sigma$)"%sig_level)

        plt.title("%d detected" %len(fsnr[idx_f]))
        plt.subplot(2,5,5)
        plt.scatter(sex_snr[idx&idx_f], fsnr[idx&idx_f], s=ms)
        plt.xlabel("SEX-SNR (%s$\sigma$)"%sig_level)
        plt.ylabel("F-SNR (%s$\sigma$)"%sig_level)

        plt.subplot(2,5,10)
        idx_s = sex_snr < 80
        plt.scatter(sex_snr[idx_s], fsnr[idx_s], s=ms)
        plt.xlabel("SEX-SNR (%s$\sigma$)"%sig_level)
        plt.ylabel("F-SNR (%s$\sigma$)"%sig_level)
    pic_nm = check_path + "check_%d.png"%rank
    if os.path.exists(pic_nm):
        os.remove(pic_nm)
    plt.savefig(pic_nm)
    plt.close()

    print(cmd, rank)
t2 = time.time()
if rank == 0:
    log = "END. operation: %s, source: %s, max radius: %.2f, stamp size: %d, filter: %s, time: %.2f, code: %s"\
          %(argv[1], total_path.split("/")[-1], max_distance, size, sex_filter, t2-t1, argv[0])
    logger.info(log)



