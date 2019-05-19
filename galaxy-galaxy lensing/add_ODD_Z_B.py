import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import time
import tool_box
import numpy
import h5py
from mpi4py import MPI


# add Z_B_MIN Z_B_MAX ODDS to the CFHT catalog

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cata_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/cfht_cata/"

t1 = time.time()

if rank == 0:
    cata_name_src = []
    names_src = os.listdir(cata_path + "field_dat/")
    for nms in names_src:
        if "w" in nms and ".dat" in nms:
            cata_name_src.append(nms.split(".")[0])
    name_src_list = tool_box.allot(cata_name_src, cpus)
    print(rank, "Before, I got %d files in %d sub-lists"%(len(cata_name_src),len(name_src_list)))

else:
    name_src_list = None

sub_src_list = comm.scatter(name_src_list, root=0)

print(rank, "I got %d files"%len(sub_src_list))

datas = []
nums = []
for i in range(1, 5):
    f = h5py.File(cata_path + "CFHT_W%d.hdf5"%i,"r")
    data = f["/data"].value
    nums.append(data.shape[0])
    datas.append(data)
    f.close()

file_header = "pos  Flag   FLUX_RADIUS  e1  e2  weight  fitclass    SNratio MASK    Z_B m   c2  LP_Mi   star_flag   MAG_i   " \
              "Z_B_MIN  Z_B_MAX ODDS"
for nms in sub_src_list:
    src_path = cata_path + "field_dat/" + nms + ".dat"
    dst_path = cata_path + "field_dat/" + nms + "_new.dat"

    # read source data
    src_data = numpy.loadtxt(src_path)
    src_sp = src_data.shape

    mask = numpy.zeros((src_sp[0], 1),dtype=numpy.intc)
    # 3 extra cols for Z_B_MIN Z_B_MAX ODDS
    dst_data = numpy.zeros((src_sp[0], src_sp[1] + 3))

    dst_data[:, :src_sp[1]] = src_data

    ra_min, ra_max = src_data[:,0].min(),src_data[:,0].max()
    dec_min, dec_max = src_data[:,1].min(),src_data[:,1].max()

    area_label = int(nms[1])
    idx1 = datas[area_label-1][:,0] >= ra_min
    idx2 = datas[area_label-1][:,0] <= ra_max
    idx3 = datas[area_label-1][:,1] >= dec_min
    idx4 = datas[area_label-1][:,1] <= dec_max

    idx = idx1 & idx2 & idx3 & idx4

    for i in range(src_sp[0]):

        ra_src, dec_src, e1_src, e2_src = src_data[i,0],src_data[i,1],src_data[i,4],src_data[i,5]

        for j in range(nums[area_label-1]):

            if idx[j]:

                ra_, dec_ = datas[area_label-1][j, 0], datas[area_label-1][j, 1]
                e1_, e2_ = datas[area_label-1][j, 2], datas[area_label-1][j, 3]

                d_ra, d_dec = numpy.abs(ra_src - ra_), numpy.abs(dec_src - dec_)
                d_e1, d_e2 = numpy.abs(e1_src - e1_), numpy.abs(e2_src - e2_)

                if d_ra < 0.00001 and d_dec < 0.00001 and d_e1 < 0.00001 and d_e2 < 0.00001:

                    dst_data[:, src_sp[1]] = datas[area_label-1][j, 5]
                    dst_data[:, src_sp[1] + 1] = datas[area_label - 1][j, 6]
                    dst_data[:, src_sp[1] + 2] = datas[area_label - 1][j, 7]

                    mask[i] = +1
    idxm = mask == 1
    if idxm.sum() != src_sp[0]:
        print(nms, "Some sources are missing")
    idxm = mask > 1
    if idxm.sum() != 0:
        print(nms, "Overlap")

    numpy.savetxt(fname=dst_path, X=dst_data, header=file_header)
t2 = time.time()
print("%.2f sec"%(t2-t1))