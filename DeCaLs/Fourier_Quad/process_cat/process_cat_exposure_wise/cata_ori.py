import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
from plot_tool import Image_Plot
from mpi4py import MPI
import tool_box
import warnings

warnings.filterwarnings('error')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


ori_cat_chara = "_all.cat"

data_band = ["r", "z"]

total_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/DECALS/cat_hdf5_ori"
ori_cat_path = "/lustre/home/acct-phyzj/share/DECALS_all_with_flat"

# if rank == 0:
#     for ib in data_band:
#         files = os.listdir(ori_cat_path + "/" + ib)
#
#         band_files = []
#
#         for fnm in files:
#             if ori_cat_chara in fnm:
#                 band_files.append("%s\n"%fnm.split(".")[0])
#
#         with open(total_path + "/code/exposures_%s_band.dat"%ib,"w") as f:
#             f.writelines(band_files)
#
# comm.Barrier()

need_idx = [0,1, 5,6, 14,21, 24,25, 26,27,28,29,30]
for ib in data_band:
    with open(total_path + "/code/exposures_%s_band.dat"%ib, "r") as f:
        cc = f.readlines()

    files = tool_box.alloc(cc, numprocs, "seq")[rank]

    exposures_candidates_avail_sub = []

    for fnm in files:
        nm = fnm.split("\n")[0]
        src_path = ori_cat_path + "/%s/%s.cat"%(ib, nm)
        dst_path = total_path + "/%s/%s.hdf5"%(ib, nm)

        try:
            src_data = numpy.loadtxt(src_path, dtype=numpy.float32)
            row, col = src_data.shape
        except:
            if os.path.exists(src_path):
                log_inform = "%d Failed in reading %s %d Bytes !\n" % (
                    rank, src_path, os.path.getsize(src_path))
                print(log_inform)
            row, col = 0, 0
        log_inform = "%d Failed in reading %s %d Bytes !\n" % (
            rank, src_path, os.path.getsize(src_path))
        print(log_inform)

        if row > 0:

            idx1 = src_data[:,15] < 48
            idx2 = src_data[:,16] < 48
            idx = idx1 & idx2
            src_num = idx.sum()
            if src_num > 0:
                dst_data = numpy.zeros((src_num, len(need_idx)), dtype=numpy.float32)
                for ii, ir in enumerate(need_idx):
                    dst_data[:,ii] = src_data[:,ir][idx]

                h5f_expo = h5py.File(dst_path, "w")
                h5f_expo["/data"] = dst_data
                h5f_expo.close()

                exposures_candidates_avail_sub.append(dst_path +"\n")

    exposures_candidates_avail = comm.gather(exposures_candidates_avail_sub, root=0)

    comm.Barrier()

    if rank == 0:

        exposures_avail_all = []
        for fsb in exposures_candidates_avail:
            exposures_avail_all.extend(fsb)

        with open(total_path + "/code/exposures_avail_%s_band.dat"%ib, "w") as f:
            f.writelines(exposures_avail_all)
    comm.Barrier()
comm.Barrier()