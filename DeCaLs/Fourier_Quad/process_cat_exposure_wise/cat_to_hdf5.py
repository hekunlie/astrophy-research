import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
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

bands = ["g", "r", "z"]

if mode == "cata_name":
    if rank == 0:

        all_band_files = []

        for ib in bands:

            band_file_path = total_path + "/%s"%id
            files_nm = os.listdir(band_file_path)

            band_files = []
            for fnm in files_nm:
                if ".cat" in fnm:
                    band_files.append(band_file_path + "/%s\n"%fnm)

            all_band_files.extend(band_files)
            print("%s band: %d exposures"%(ib, len(band_files)))
            with open("nname_exposures_%s.dat"%ib, "w") as f:
                f.writelines(band_files)

        with open("nname_exposures_all.dat", "w") as f:
            f.writelines(all_band_files)
        print("All %d exposures"%len(all_band_files))

elif mode == "hdf5_cata":
    # convert the .dat to .hdf5

    all_exposures = []

    with open("nname_exposures_all.dat", "r") as f:
        contents = f.readlines()
    for cc in contents:
        all_exposures.append(cc.split(".")[0])

    sub_exposures = tool_box.alloc(all_exposures, cpus, "seq")[rank]

    src_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/DECALS"

    sub_exposures_avail = []
    exception_sub = []
    for expo_nm in sub_exposures:
        # read the the field data
        exposure_name = expo_nm.split("/")[-1]
        band_name = expo_nm.split("/")[-2]

        expo_src_path = expo_nm + ".cat"
        expo_h5_path = src_path + "/%s/%s.hdf5"%(band_name, exposure_name)

        try:
            edat = numpy.loadtxt(expo_src_path, dtype=numpy.float32)

            # Nan check
            idx = numpy.isnan(edat)
            if idx.sum() > 0:
                log_inform = "Find NAN in %s"%expo_src_path
                exception_sub.append(log_inform)
                print(log_inform)

            h5f_expo = h5py.File(expo_h5_path, "w")
            h5f_expo["/data"] = edat
            h5f_expo.close()
            sub_exposures_avail.append(expo_h5_path + "\n")

        except:
            if os.path.exists(expo_src_path):
                log_inform = "%d Failed in reading %s %d Bytes !\n" % (
                rank, expo_src_path, os.path.getsize(expo_src_path))
            else:
                log_inform = "%d can't find %s!\n" % (rank, expo_src_path)
            exception_sub.append(log_inform)

    exposures_collection = comm.gather(sub_exposures_avail, root=0)
    exception_collection = comm.gather(exception_sub, root=0)

    comm.Barrier()
    if rank == 0:
        exception_all = []
        for ec in exception_collection:
            exception_all.extend(ec)

        exposures_avail = []
        for fsb in exposures_collection:
            exposures_avail.extend(fsb)
        with open("nname_exposures_all_avail.dat","w") as f:
            f.writelines(exposures_avail)

        exception_all.append("Totally: %d exposures\n"%(len(exposures_avail)))
        with open("log.dat","w") as f:
            f.writelines(exception_all)