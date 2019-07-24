import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import time
import tool_box
import numpy
import h5py
from mpi4py import MPI


#############################################################################################
# add Z_B_MIN Z_B_MAX ODDS (from the .csv files) to the CFHT catalog
# it will create two new files (.hdf5 & _new.dat) that contains the added parameters
# the hdf5 file contains 3 arrays:
# 1. "data": the catalog with the 3 parameters
# 2. "mask": it should be 1 for each source
# 3. "dRA_dDEC": delta RA and delta DEC, they should be very small for each source ( < 1e-5)
############################################################################################


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

# cata_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/cfht_cata/"
# log_path = "/home/hkli/work/test/log/log_%d.dat"%rank

cata_path = "/mnt/perc/hklee/CFHT/catalog/cfht_cata/"
log_path = "/home/hklee/work/test/log/log_%d.dat"%rank

logger = tool_box.get_logger(log_path)

t1 = time.time()

if rank == 0:
    cata_name_src = []
    names_src = os.listdir(cata_path + "field_dat/")
    for nms in names_src:
        if "w" in nms and ".dat" in nms and "_new" not in nms and "-pz" not in nms:
            cata_name_src.append(nms.split(".")[0])
    name_src_list = tool_box.allot(cata_name_src, cpus)
    print(rank, "Before, I got %d files in %d sub-lists"%(len(cata_name_src),len(name_src_list)))

else:
    name_src_list = None

sub_src_list = comm.scatter(name_src_list, root=0)

print(rank, "I got %d files"%len(sub_src_list))

# data is the catalog contains additional parameters, "###.tsv"
# read all data from each area
for i in range(1, 5):
    h5f = h5py.File(cata_path + "CFHT_W%d.hdf5"%i,"r")
    temp = h5f["/data"].value
    h5f.close()
    if i == 1:
        data = temp
    else:
        data = numpy.row_stack((data, temp))

num = data.shape[0]

file_header = "RA   DEC Flag   FLUX_RADIUS  e1  e2  weight  fitclass    SNratio    MASK   Z_B " \
              "  m   c2  LP_Mi   star_flag   MAG_i   Z_B_MIN  Z_B_MAX  ODDS"

for nms in sub_src_list:

    st1 = time.time()

    logger.info("Start %s"%nms)
    src_path = cata_path + "field_dat/" + nms + ".dat"
    dst_path = cata_path + "field_dat/" + nms + "_new.dat"
    dst_path_h5 = cata_path + "field_dat/" + nms + ".hdf5"

    # read source data
    src_data = numpy.loadtxt(src_path)
    src_sp = src_data.shape

    mask = numpy.zeros((src_sp[0], ),dtype=numpy.intc) - 99
    # delta RA & delta DEC
    delta_position = numpy.zeros((src_sp[0], 2)) - 99
    # 3 extra cols for Z_B_MIN Z_B_MAX ODDS
    dst_data = numpy.zeros((src_sp[0], src_sp[1] + 3), dtype=numpy.float64) - 99.

    dst_data[:, :src_sp[1]] = src_data

    ra_min, ra_max = src_data[:,0].min(),src_data[:,0].max()
    dec_min, dec_max = src_data[:,1].min(),src_data[:,1].max()

    area_label = int(nms[1])
    idx1 = data[:,0] >= ra_min - 0.001
    idx2 = data[:,0] <= ra_max + 0.001
    idx3 = data[:,1] >= dec_min - 0.001
    idx4 = data[:,1] <= dec_max + 0.001

    idx = idx1 & idx2 & idx3 & idx4

    logger.info("Begin to match...")
    sub_data = data[idx]
    # loop the source in original catalog to find the target in the new download catalog
    for i in range(src_sp[0]):

        ra_src, dec_src = src_data[i,0], src_data[i,1]
        # radius
        del_radius = numpy.abs(ra_src - sub_data[:,0]) + numpy.abs(dec_src - sub_data[:,1])

        del_radius_min = del_radius.min()

        if del_radius_min <= 0.00001:

            npw_dr = numpy.where(del_radius == del_radius_min)[0][0]

            dst_data[i, src_sp[1]] = sub_data[npw_dr, 5]
            dst_data[i, src_sp[1] + 1] = sub_data[npw_dr, 6]
            dst_data[i, src_sp[1] + 2] = sub_data[npw_dr, 7]
            # delta ra & delta dec for checking
            delta_position[i, 0] = sub_data[npw_dr, 0] - ra_src
            delta_position[i, 1] = sub_data[npw_dr, 1] - dec_src
            mask[i] = +1

    st2 = time.time()

    idxm = mask == 1
    if idxm.sum() != src_sp[0]:
        print(nms, "Some sources are missing. source: %d. matched: %d. diff: %d"%(src_sp[0],idxm.sum(), src_sp[0] - idxm.sum()))
        logger.info("Done %s,Some sources are missing. %.2f sec" % (nms, st2 - st1))
    idxm = mask > 1
    if idxm.sum() != 0:
        print(nms, "Overlap")
        logger.info("Done %s,Some Overlap. %.2f sec" % (nms, st2 - st1))

    numpy.savetxt(fname=dst_path, X=dst_data, header=file_header)
    h5f = h5py.File(dst_path_h5, "w")
    h5f["/data"] = dst_data
    h5f["/dRA_dDEC"] = delta_position
    h5f["/mask"] = mask
    h5f.close()

    logger.info("write file %s. MASK: max = %d, min = %d. Delta R: max = %.6f"%(dst_path, mask.max(), mask.min(), delta_position.max()))

t2 = time.time()
print("%.2f sec"%(t2-t1))