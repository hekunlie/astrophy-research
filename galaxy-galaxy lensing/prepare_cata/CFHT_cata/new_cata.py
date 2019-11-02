import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
import time
import tool_box
import numpy
import h5py
from mpi4py import MPI


#############################################################################################
# segment the newly-download CFHT catalog according to the current catalog.
# the hdf5 file contains 3 arrays:
# 1. "data": the catalog with the 3 parameters
# 2. "mask": it should be 1 for each source
# 3. "dRA_dDEC": delta RA and delta DEC, they should be very small for each source ( < 1e-5)
############################################################################################


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


cata_path = "/mnt/perc/hklee/CFHT/catalog/cfht_cata/"
log_path = "/home/hklee/work/test/log/log_%d.dat"%rank

logger = tool_box.get_logger(log_path)

t1 = time.time()

for area_id in range(int(argv[1]),int(argv[1])+1):

    if rank == 0:

        cata_name_src = []

        # read the whole field names
        with open(cata_path+"nname.dat","r") as f:
            field_names = f.readlines()

        for nm in field_names:
            if "w%d"%area_id in nm:
                cata_name_src.append(nm.split("\n")[0])

        name_src_list = tool_box.allot(cata_name_src, cpus)
        print("Area: %d has %d fields"%(area_id,len(cata_name_src)))

    else:
        name_src_list = None

    sub_src_list = comm.scatter(name_src_list, root=0)

    old_file_header = "RA   DEC  Flag   FLUX_RADIUS  e1  e2  weight  fitclass   " \
                      "SNratio  MASK Z_B m   c2  LP_Mi   star_flag   MAG_i"

    new_file_header = "RA   DEC     Level   Flag    FLUX_RADIUS     CLASS_STAR      e1      e2      " \
                      "weight  fitclass    SNratio MASK    Z_B     Z_B_MIN Z_B_MAX ODDS    m       c2      LP_Mi   " \
                      "star_flag       MAG_u   MAG_g   MAG_r   MAG_i   MAG_y   MAG_z"

    src_index = old_file_header.split()
    dst_index = new_file_header.split()

    src_ra_id, src_dec_id = src_index.index("RA"), src_index.index("DEC")
    dst_ra_id, dst_dec_id = dst_index.index("RA"), dst_index.index("DEC")

    src_zb_id, dst_zb_id = src_index.index("Z_B"), dst_index.index("Z_B")

    if rank == 0:
        print(src_ra_id, src_dec_id,dst_ra_id, dst_dec_id, src_zb_id, dst_zb_id)


    h5f = h5py.File(cata_path + "new/w_%d.hdf5"%area_id, "r")
    area_data = h5f["/data"].value
    h5f.close()

    comm.Barrier()

    for nms in sub_src_list:

        st1 = time.time()

        logger.info("Start %s"%nms)
        src_path = cata_path + "field_dat/" + nms + ".dat"
        dst_path = cata_path + "field_dat_new/" + nms + ".dat"
        dst_path_h5 = cata_path + "field_dat_new/" + nms + ".hdf5"

        # read old catalog
        src_data = numpy.loadtxt(src_path)
        src_sp = src_data.shape

        mask = numpy.zeros((src_sp[0], ), dtype=numpy.intc)
        # delta RA, delta DEC, delta Z_B
        diff_check = numpy.zeros((src_sp[0], 3))

        # for new catalog
        dst_row = src_sp[0]
        dst_cols = len(dst_index)
        dst_data = numpy.zeros((dst_row,dst_cols))

        ra_min, ra_max = src_data[:,src_ra_id].min(),src_data[:,src_ra_id].max()
        dec_min, dec_max = src_data[:,src_dec_id].min(),src_data[:,src_dec_id].max()

        idx1 = area_data[:,dst_ra_id] >= ra_min - 0.005
        idx2 = area_data[:,dst_ra_id] <= ra_max + 0.005
        idx3 = area_data[:,dst_dec_id] >= dec_min - 0.005
        idx4 = area_data[:,dst_dec_id] <= dec_max + 0.005

        idx = idx1 & idx2 & idx3 & idx4

        logger.info("Begin to match...")

        sub_area_data = area_data[idx]
        # loop the source in old catalog to find the target in the new download catalog
        for i in range(dst_row):
            ra_src, dec_src = src_data[i,src_ra_id], src_data[i,src_dec_id]

            # difference
            del_radius = numpy.abs(ra_src - sub_area_data[:,dst_ra_id]) + numpy.abs(dec_src - sub_area_data[:,dst_dec_id])

            del_radius_min = del_radius.min()

            if del_radius_min <= 0.00001:

                npw_dr = numpy.where(del_radius == del_radius_min)[0][0]

                dst_data[i] = sub_area_data[npw_dr]

                # delta ra & delta dec for checking
                diff_check[i, 0] = numpy.abs(sub_area_data[npw_dr, dst_ra_id] - ra_src)
                diff_check[i, 1] = numpy.abs(sub_area_data[npw_dr, dst_dec_id] - dec_src)
                diff_check[i, 2] = numpy.abs(sub_area_data[npw_dr, dst_zb_id] - src_data[i,src_zb_id])
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

        numpy.savetxt(fname=dst_path, X=dst_data, header=new_file_header)

        h5f = h5py.File(dst_path_h5, "w")
        h5f["/data"] = dst_data
        h5f["/diff"] = diff_check
        h5f["/mask"] = mask
        h5f.close()

        logger.info("write file %s. MASK: max = %d, min = %d. Delta R: max = %.6f"%(dst_path, mask.max(), mask.min(), diff_check.max()))

    comm.Barrier()
    t2 = time.time()
    print("%.2f sec"%(t2-t1))