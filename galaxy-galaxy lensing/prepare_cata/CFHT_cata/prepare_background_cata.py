import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
import h5py
from mpi4py import MPI
from sys import argv
import numpy
import time
from subprocess import Popen

################################################################################################
# collect: collect the data from the files of each field. It creates the "cfht_cata.hdf5" in
#           the parent directory of the one contain the field catalog.
#           If the catalog file doesn't exist, run it firstly !!!.
#           CFHT catalog contains 19 (0~18) columns,  19'th & 20'th column are the PZ data from Dong FY.
#
# select: select the galaxy
#
################################################################################################

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]
cmds = ["collect", "select"]

if cmd not in cmds:
    if rank == 0:
        print("parameter must be one of ", cmds)
    exit()

area_num = 4

cata_path = "/mnt/perc/hklee/CFHT/catalog/"
data_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/"

cfht_cata_path = cata_path + "cfht_cata/"


############################# CFHTLenS Option ####################################################
# CFHT catalog contains 19 (0~18) columns,  19'th & 20'th column are the PZ data from Dong FY.
ra_lb_c = 0
dec_lb_c = 1
flag_lb_c = 2
flux_rad_lb_c = 3
e1_lb_c = 4
e2_lb_c = 5
weight_lb_c = 6
fitclass_lb_c = 7
snr_lb_c = 8
mask_lb_c = 9
z_lb_c = 10
m_lb_c = 11
c_lb_c = 12
lp_mi_lb_c = 13
starflag_lb_c = 14
mag_lb_c = 15
z_min_lb_c = 16
z_max_lb_c = 17
odds_lb_c = 18
pz1_lb = 19
pz2_lb = 20

# change the cutoff threshold in the part of "select_cfht"

############################# CFHTLenS Option ####################################################

logger = tool_box.get_logger("%s/work/test/log/%d.dat"%(my_home, rank))


############################# CFHTLenS catalog collection ########################################
# combine the data of each field into one big catalog
if cmd == "collect":

    t1 = time.time()

    logger.info("RANK: %03d. Begin..."%rank)

    h5f_path = cfht_cata_path + "cfht_cata.hdf5"
    if rank == 0:
        h5f = h5py.File(h5f_path, "w")
        h5f.close()

    pre_fields = tool_box.field_dict(cfht_cata_path + "nname.dat")[1]

    # loop the areas,
    # the catalog of each area will be stored in "w_i"
    num = []
    for area_id in range(1,area_num+1):
        # distribution the files
        field_tar = []
        for field in pre_fields:
            if "w%d" % area_id in field:
                field_tar.append(field)

        # if fields are more than the threads,
        # some thread will get an empty "field_pool"
        field_pool = tool_box.allot(field_tar, cpus)[rank]

        # check
        anothers = ["w%d"%i for i in range(1,area_num+1) if i != area_id]
        for field_name in field_pool:
            for ext_field in anothers:
                if ext_field in field_name:
                    print("WRONG",rank, area_id, field_name)
                    exit()

        field_count = 0

        for field_name in field_pool:
            # c_data_path = cfht_cata_path + "field_dat/%s_new.dat"%field_name
            c_data_path = cfht_cata_path + "field_dat/%s.hdf5"%field_name
            pz_data_path = cfht_cata_path + "field_dat/%s-pz.dat"%field_name

            if os.path.exists(c_data_path):
                # try:

                # c_data = numpy.loadtxt(c_data_path)
                h5f_c = h5py.File(c_data_path,"r")
                c_data = h5f_c["/data"].value
                h5f_c.close()

                pz_data = numpy.loadtxt(pz_data_path)

                point_data = numpy.column_stack((c_data, pz_data))

                if field_count == 0:
                    cata_data = point_data
                else:
                    cata_data = numpy.row_stack((cata_data, point_data))
                field_count += 1

                # except:
                #     print(rank, "%s.dat doesn't exist"%field_name)
            else:
                print("Rank %d Can't find %s"%(rank,field_name))

        # in case of that some thread gets nothing from the file (non-existed)
        if field_count < 1:
            cata_data = numpy.zeros((1,1))
        data_sp = cata_data.shape
        data_sps = comm.gather(data_sp, root=0)
        num.append(data_sp)


        # DEBUG
        # for ir in range(cpus):
        #     if ir == rank:
        #         print(rank, len(field_tar), field_pool, field_count, cata_data.shape, type(cata_data), cata_data.dtype)
        #         h5f_temp = h5py.File("%d.hdf5"%rank,"w")
        #         h5f_temp["/data"] = cata_data
        #         h5f_temp.close()
        #     comm.Barrier()

        if rank > 0:
            logger.info("RANK: %03d. AREA: %d. Send... Files: %d" % (rank, area_id, field_count))
            comm.Send([cata_data,MPI.DOUBLE], dest=0, tag=rank)
        else:
            logger.info("RANK: %03d. AREA: %d. Receive... Files: %d" % (rank, area_id, field_count))
            if data_sp[0] > 1 and data_sp[1] > 1:
                stack_pool = [cata_data]
            else:
                stack_pool = []

            for procs in range(1,cpus):
                recvs = numpy.empty(data_sps[procs], dtype=numpy.float64)
                comm.Recv(recvs, source=procs, tag=procs)
                if data_sps[procs][0] > 1 and data_sps[procs][1] > 1:
                    stack_pool.append(recvs)

            for stack_count, arr in enumerate(stack_pool):
                if stack_count == 0:
                    recv_buffer = arr
                else:
                    recv_buffer = numpy.row_stack((recv_buffer, arr))

            h5f = h5py.File(h5f_path,"r+")

            h5f["/w_%d" % area_id] = recv_buffer
            h5f.close()
            logger.info("RANK: %03d. AREA: %d. RECV: (%d, %d)" % (rank, area_id, recv_buffer.shape[0],recv_buffer.shape[1]))
            print("RECV: ", recv_buffer.shape)
        comm.Barrier()

    print(rank, "the galaxy number in each area", num)
    log_str = "RANK: %03d. galaxy number in each area [" % rank
    for i in range(0, area_num):
        log_str += "(%d, %d), "%(num[i][0], num[i][1])
    log_str += "]. "
    logger.info(log_str)

    # stack the sub-catalogs from each area
    t2 = time.time()
    if rank == 0:

        h5f = h5py.File(h5f_path,"r+")
        for area_id in range(1,area_num+1):
            logger.info("RANK: %03d. Stacking %d..." % (rank, area_id))
            temp_s = h5f["/w_%d" % area_id].value
            if area_id == 1:
                data = temp_s
            else:
                data = numpy.row_stack((data, temp_s))
            print("Totally, %d galaxies are detected in W_%d" % (len(temp_s), area_id))

        h5f["/total"] = data
        h5f.close()
        logger.info("RANK: %03d. Finish..." %rank)
    t3 = time.time()
    if rank == 0:
        print("%.2f,%.2f, %d"%(t2-t1,t3-t2, len(data)))
############################# CFHTLenS catalog collection ########################################


############################# CFHTLenS catlog cutoff #############################################
if cmd == "select":
    t1 = time.time()

    block_scale = float(argv[2])
    margin = 0.1*block_scale

    h5f_path = cfht_cata_path + "cfht_cata.hdf5"
    h5f_path_cut = data_path + "cfht_cata/cfht_cata_cut.hdf5"

    if rank == 0:
        h5f = h5py.File(h5f_path_cut, "w")
        h5f["/block_scale"] = numpy.array([block_scale],dtype=numpy.float64)
        h5f.close()
    comm.Barrier()

    if rank < area_num:
        h5f = h5py.File(h5f_path, "r")
        cata_data = h5f["/w_%d"%(rank+1)].value
        h5f.close()

        # cut off
        epsilon = 0.0000001

        # weight > 0
        idx_weight = cata_data[:, weight_lb_c] > 0 - epsilon

        # mask <= 1
        idx_mask = cata_data[:, mask_lb_c] <= 1 + epsilon

        # FITCLASS = 0
        idx_fitclass = numpy.abs(cata_data[:, fitclass_lb_c]) < 0 + epsilon

        # MAG < 24.7
        idx_mag = cata_data[:, mag_lb_c] < 24.7

        # Redshift
        idxz_1 = cata_data[:, z_lb_c] <= 15
        idxz_2 = cata_data[:, z_lb_c] >= 0
        idxz_3 = cata_data[:, odds_lb_c] > 0.5 - epsilon

        # # the PZ selection criteria from Dong FY
        # idx_pz1 = numpy.abs(cata_data[:, pz1_lb] - 1) <= epsilon
        # idx_pz2 = cata_data[:, pz2_lb] <= 0.2

        cut_idx = idx_weight & idx_mask & idx_fitclass & idx_mag & idxz_1 & idxz_2 & idxz_3

        ra = cata_data[:, ra_lb_c][cut_idx]
        dec = cata_data[:, dec_lb_c][cut_idx]
        cos_dec = numpy.abs(numpy.cos((dec/180.*numpy.pi)))

        redshift = cata_data[:, z_lb_c][cut_idx]
        z_min = cata_data[:,z_min_lb_c][cut_idx]
        z_max = cata_data[:,z_max_lb_c][cut_idx]
        odds = cata_data[:,odds_lb_c][cut_idx]

        mag = cata_data[:, mag_lb_c][cut_idx]

        m_bias = cata_data[:,m_lb_c][cut_idx]
        c_bias = cata_data[:,c_lb_c][cut_idx]

        e1 = cata_data[:,e1_lb_c][cut_idx]
        # the minus arises from the true RA-axis is in opposite direction
        e2 = -(cata_data[:,e2_lb_c][cut_idx] - c_bias)

        weight = cata_data[:,weight_lb_c][cut_idx]

        # starflag = cata_data[:, starflag_lb_c][cut_idx]
        #
        # fitclass = cata_data[:,fitclass_lb_c][cut_idx]
        #
        # mask = cata_data[:,mask_lb_c][cut_idx]

        data_num = len(redshift)

        names = ["Z", "RA", "DEC", "MAG", "COS_DEC", "Z_MIN", "Z_MAX", "ODDS",
                 "E1", "E2", "WEIGHT", "M", "C"]

        datas = [redshift, ra, dec, mag, cos_dec, z_min, z_max, odds,
                 e1, e2, weight, m_bias, c_bias]

        # set up RA & DEC bin
        ra_min, ra_max = ra.min()-margin, ra.max()+margin
        dec_min, dec_max = dec.min()-margin, dec.max()+margin

        nx = int((ra_max-ra_min)/block_scale) + 2
        ny = int((dec_max-dec_min)/block_scale) + 2
        grid_num = nx*ny
        grid_shape = numpy.array([ny, nx], dtype=numpy.intc)

        ra_bin = numpy.zeros((nx+1, 1))
        dec_bin = numpy.zeros((ny+1, 1))
        for i in range(nx+1):
            ra_bin[i] = ra_min + i*block_scale
        for i in range(ny+1):
            dec_bin[i] = dec_min + i*block_scale
        if ra_bin.max() < ra_max:
            print("Too less RA bins")
            exit(0)
        if dec_bin.max() < dec_max:
            print("Too less DEC bins")
            exit(0)

        # the boundary of each block
        boundx = numpy.zeros((grid_num, 4))
        boundy = numpy.zeros((grid_num, 4))
        for i in range(ny):
            ix = i * nx
            for j in range(nx):
                tag = ix + j
                boundy[tag, 0] = dec_bin[i]
                boundy[tag, 1] = dec_bin[i]
                boundy[tag, 2] = dec_bin[i + 1]
                boundy[tag, 3] = dec_bin[i + 1]

                boundx[tag, 0] = ra_bin[j]
                boundx[tag, 1] = ra_bin[j + 1]
                boundx[tag, 2] = ra_bin[j]
                boundx[tag, 3] = ra_bin[j + 1]

        # # galaxy count in each block
        # num_in_block = numpy.zeros((1, grid_num), dtype=numpy.intc)
        # block_start = numpy.zeros((1, grid_num), dtype=numpy.intc)
        # # the galaxy labels in the block
        # gal_sequence = numpy.zeros((1, data_num), dtype=numpy.intc)

        # galaxy count in each block
        num_in_block = numpy.zeros((grid_num,), dtype=numpy.intc)
        block_start = numpy.zeros((grid_num,), dtype=numpy.intc)
        # the galaxy labels in the block
        gal_in_block = numpy.zeros((data_num,), dtype=numpy.intc)
        # block label of each gal
        block_label = numpy.zeros((data_num, ), dtype=numpy.intc)
        # the galaxy label
        gal_label = numpy.arange(0, data_num)

        for i in range(ny):
            idx_1 = dec >= dec_bin[i]
            idx_2 = dec < dec_bin[i+1]
            idx_sub = idx_1 & idx_2
            ix = i*nx
            for j in range(nx):
                tag = ix + j
                idx_3 = ra >= ra_bin[j]
                idx_4 = ra < ra_bin[j+1]
                idx = idx_3 & idx_4 & idx_sub
                num_in_block[tag] = idx.sum()

                block_start[tag] = num_in_block[:tag].sum()

                gal_in_block[block_start[tag]: block_start[tag]+num_in_block[tag]] = gal_label[idx]
                block_label[idx] = tag

        print("W%d: %d ==> %d"%(rank, cata_data.shape[0], data_num))
    comm.Barrier()

    for area_id in range(area_num):
        if rank == area_id:
            h5f = h5py.File(h5f_path_cut, "r+")

            h5f["/w_%d/grid_shape" % (rank + 1)] = grid_shape

            h5f["/w_%d/num_in_block" % (rank + 1)] = num_in_block
            h5f["/w_%d/block_start" % (rank + 1)] = block_start
            h5f["/w_%d/gal_in_block" % (rank + 1)] = gal_in_block
            h5f["/w_%d/block_label" % (rank + 1)] = block_label

            h5f["/w_%d/block_boundy" % (rank + 1)] = boundy
            h5f["/w_%d/block_boundx" % (rank + 1)] = boundx
            h5f["/w_%d/RA_bin" % (rank + 1)] = ra_bin
            h5f["/w_%d/DEC_bin" % (rank + 1)] = dec_bin

            for i in range(len(names)):
                h5f["/w_%d/%s"%(rank+1,names[i])] = datas[i]
            h5f.close()
        comm.Barrier()
    t2 = time.time()
    if rank == 0:
        for i in range(area_num):
            cmd = "../add_com_dist %s /w_%d/" % (h5f_path_cut, i + 1)
            a = Popen(cmd, shell=True)
            a.wait()
        print("%s, %.2f sec"%(tool_box.get_time_now(),t2-t1))
############################# CFHTLenS catlog cutoff #############################################