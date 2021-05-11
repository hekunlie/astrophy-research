import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
import warnings

warnings.filterwarnings('error')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

total_path = argv[1]
mode = argv[2]
ori_cat_path = argv[3]

ori_cat_chara = "_all_alt.cat"

data_band = ["r", "z"]

if mode == "cata_name":
    if rank == 0:

        # with open("anomaly_expo.dat", "r") as f:
        #     contents = f.readlines()
        # anomaly_expo = []
        # anomaly_val = []
        # for cc in contents:
        #     val, field_nm, expo = cc.rstrip().split()[1:4]
        #     anomaly_expo.append(int(expo.split("p")[0]))
        #     anomaly_val.append(float(val))
        #     # print(anomaly_expo[-1],anomaly_val[-1])
        all_files = []
        for band in data_band:
            files_nm = os.listdir(ori_cat_path + "/%s"%(band))
            band_files = []

            for fnm in files_nm:
                if ori_cat_chara in fnm:
                    band_files.append("%s\n"%fnm.split(".")[0])
                    all_files.append("%s %s\n"%(fnm.split(".")[0],band))

            with open(total_path + "/cat_inform/exposure_name_%s_band.dat"%band, "w") as f:
                f.writelines(band_files)

            print(len(band_files), " exposures in %s band\n%s/%s"%(band,ori_cat_path,band))

        with open(total_path + "/cat_inform/exposure_name_all_band.dat", "w") as f:
            f.writelines(all_files)
        print(len(all_files), " exposures")

if mode == "hdf5_cata":
    # convert the .dat to .hdf5

    exposures_candidates = []
    exposures_candidates_band = []
    with open(total_path + "/cat_inform/exposure_name_all_band.dat", "r") as f:
        contents = f.readlines()
    for expo_name in contents:
        expo_nm, band = expo_name.split("\n")[0].split()
        exposures_candidates.append(expo_nm)
        exposures_candidates_band.append(band)

    exposures_candidates_sub = tool_box.alloc(exposures_candidates, cpus, "seq")[rank]
    exposures_candidates_band_sub = tool_box.alloc(exposures_candidates_band, cpus, "seq")[rank]

    exposures_candidates_avail_sub = []
    exception_sub = []

    for tag, fns in enumerate(exposures_candidates_sub):
        # read the the field data
        iband = exposures_candidates_band_sub[tag]
        expo_src_path = ori_cat_path + "/%s/%s.cat" % (iband,fns)
        expo_h5_path = total_path + "/cat_hdf5/%s/%s.hdf5" % (iband,fns)
        try:
            src_data = numpy.loadtxt(expo_src_path, dtype=numpy.float32)
            row, col = src_data.shape
        except:
            if os.path.exists(expo_src_path):
                log_inform = "%d Failed in reading %s %d Bytes !\n" % (
                    rank, expo_src_path, os.path.getsize(expo_src_path))
            else:
                log_inform = "%d can't find %s!\n" % (rank, expo_src_path)
            exception_sub.append(log_inform)
            row, col = 0, 0

        if row > 0:
            # Nan check
            idx = numpy.isnan(src_data)
            if idx.sum() > 0:
                print("Find Nan in ", expo_src_path)

            h5f_expo = h5py.File(expo_h5_path, "w")
            h5f_expo["/data"] = src_data
            h5f_expo.close()

            exposures_candidates_avail_sub.append(expo_h5_path + " %s"%exposures_candidates_band_sub[tag])

    exposures_candidates_avail = comm.gather(exposures_candidates_avail_sub, root=0)
    exception_collection = comm.gather(exception_sub, root=0)

    comm.Barrier()

    if rank == 0:
        exception_all = []
        for ec in exception_collection:
            exception_all.extend(ec)

        exposures_collect = []
        for fsb in exposures_candidates_avail:
            exposures_collect.extend(fsb)

        exposures_avail_all = []
        exposures_avail_band = [[] for i in range(len(data_band))]

        for fsb in exposures_collect:
            expo_path, band = fsb.split()
            exposures_avail_all.append(expo_path + "\n")

            for i in range(len(data_band)):
                if band == data_band[i]:
                    exposures_avail_band[i].append(expo_path + "\n")

        with open(total_path + "/cat_inform/exposure_avail_all.dat","w") as f:
            f.writelines(exposures_avail_all)

        for tag, band in enumerate(data_band):
            with open(total_path + "/cat_inform/exposure_avail_%s_band.dat"%band, "w") as f:
                f.writelines(exposures_avail_band[tag])

            exception_all.append("%d available exposures in %s band"%(len(exposures_avail_band[tag]), band))

        exception_all.append("Totally: %d available exposures\n"%len(exposures_avail_all))
        with open("log.dat", "w") as f:
            f.writelines(exception_all)
