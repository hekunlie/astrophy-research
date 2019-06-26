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


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]
cmds = ["collect_cfht", "select_cfht", "collect_fourier", "select_fourier"]

# collect: collect the data from the files of each field. run it firstly.
#           add MAG, Z, Z_MIN, Z_MAX, ODDS, starflag, weight,
#           mask, e1, e2, m, c, and fitclass to end of each row in Fourier catalog

# select: select the galaxy

if cmd not in cmds:
    if rank == 0:
        print("parameter must be one of ", cmds)
    exit()

area_num = 4

result_source = "result_ext"

envs_path = "%s/work/envs/envs.dat"%my_home

gets_item = [["cfht", "cfht_path_catalog", "0"], ["gg_lensing", "ggl_path_data", "0"]]
path_items = tool_box.config(envs_path, ["get", "get"], gets_item)

cata_path, data_path = path_items

cfht_cata_path = cata_path + "cfht_cata/"
fourier_cata_path = cata_path + "fourier_cata/"


gets_item = [["fresh_para_idx", "nstar", "0"], ["fresh_para_idx", "flux_alt", "0"],
             ["fresh_para_idx", "ra", "0"], ["fresh_para_idx", "dec", "0"],
             ["fresh_para_idx", "gf1", "0"], ["fresh_para_idx", "gf2", "0"],
             ["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
             ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
             ["fresh_para_idx", "h2", "0"], ["fresh_para_idx", "total_area", "0"]]
gets = ["get" for i in range(len(gets_item))]
para_items = tool_box.config(envs_path, gets, gets_item)

############################# Fourier_Quad Option ################################################
# the column label in Fourier_Quad catalog,
# used in Fourier_quad related option
nstar_lb = int(para_items[0])
flux_alt_lb = int(para_items[1])
total_area_lb = int(para_items[11])

ra_lb = int(para_items[2])
dec_lb = int(para_items[3])

field_g1_lb = int(para_items[4])
field_g2_lb = int(para_items[5])

mg1_lb = int(para_items[6])
mg2_lb = int(para_items[7])
mn_lb = int(para_items[8])
mu_lb = int(para_items[9])
mv_lb = int(para_items[10])

# it would change as more parameters added into the file
# the output of Fourier_quad pipeline contains 21 columns,
# but there is an additional useless columns (22'th) added before...
data_col = 21

add_col = 6

starflag_lb_f = 21
mag_lb_f = 22
z_lb_f = 23
z_min_lb_f = 24
z_max_lb_f = 25
odds_lb_f = 26
# optional
fitclass_lb_f = 27
weight_lb_f = 28
mask_lb_f = 29
e1_lb_f = 30
e2_lb_f = 31
m_lb_f = 32
c_lb_f = 33

# cut off
flux_alt_thresh = 3.25
nstar_thresh = 12
total_area_thresh = 1
field_g1_bound = 0.005
field_g2_bound = 0.0075
z_min, z_max = 0.0, 15
c1_correction = -0.000578
c2_correction = 0.000493

# flux_alt_thresh = 0
# nstar_thresh = 12
# total_area_thresh = 1
# field_g1_bound = 0.2
# field_g2_bound = 0.2
# z_min, z_max = 0.0, 15
# c1_correction = -0.000578
# c2_correction = 0.000493
############################# Fourier_Quad Option ################################################


############################# CFHTLenS Option ####################################################
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

# change the cutoff threshold in the part of "select_cfht"

############################# CFHTLenS Option ####################################################

logger = tool_box.get_logger("./logs/%d.dat"%rank)

############################# Fourier_quad data collection #######################################
# combine the data of each field into one big catalog
if cmd == "collect_fourier":

    t1 = time.time()

    h5f_path = data_path + "fourier_cata_%s.hdf5" % result_source
    if rank == 0:
        h5f = h5py.File(h5f_path, "w")
        h5f.close()

    c_dicts, c_fields = tool_box.field_dict(cfht_cata_path + "nname.dat")
    f_dicts, f_fields = tool_box.field_dict(fourier_cata_path + "nname.dat")

    # choose the fields both exist in the two catalog
    pre_fields = []
    for field in c_fields:
        if field in f_fields:
            pre_fields.append(field)

    # loop the areas,
    # the catalog of each area will be stored in "w_i"
    num = []
    for area_id in range(1, area_num+1):
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
        cata_data = None
        for field_name in field_pool:
            # read the catalogs from Fourier_Quad and CFHT
            f_data_path = fourier_cata_path + "%s/%s/%s_shear.dat"%(field_name, result_source, field_name)
            # c_data_path = cfht_cata_path + "%s.dat"%field_name
            c_data_path = cfht_cata_path + "field_dat/%s.hdf5"%field_name

            if os.path.exists(f_data_path) and os.path.exists(c_data_path):
                try:
                    # the file in Fourier_Quad catalog may be empty
                    f_data = numpy.loadtxt(f_data_path)
                    f_sp = f_data.shape
                    # the redshift magnitude and starflag, if non-exist in the cfht catalog, labeled by -1
                    new_data = numpy.zeros((f_sp[0], data_col + add_col), dtype=numpy.double) - 1
                    # the output of Fourier_quad pipeline contains 21 columns,
                    # but there is an additional useless columns (22'th) added before...
                    new_data[:,:data_col] = f_data[:,:data_col]

                    # c_data = numpy.loadtxt(c_data_path)
                    h5f_c = h5py.File(c_data_path,"r")
                    c_data = h5f_c["/data"].value
                    h5f_c.close()

                    c_sp = c_data.shape
                    # f_sp[0] > c_sp[0]
                    pairs = 0
                    mask = numpy.zeros((f_sp[0],), dtype=numpy.intc)
                    for i in range(f_sp[0]):
                        # fortran start from 1 not 0
                        # the row labels corresponding to the galaxies in cfht catalog
                        # every source in Fourier_Quad ext_cata has label corresponding
                        # to the row in the CFHTLenS catalog
                        tag = int(f_data[i, 0] - 1)
                        # check
                        del_radius = numpy.abs(f_data[i,ra_lb] - c_data[tag, 0]) + numpy.abs(f_data[i,dec_lb] - c_data[tag, 1])

                        if del_radius < 0.00005:
                            # starflag
                            new_data[i, starflag_lb_f] = c_data[tag, 14]
                            # magnitude
                            new_data[i, mag_lb_f] = c_data[tag, 15]
                            # redshift
                            new_data[i, z_lb_f] = c_data[tag, 10]
                            # Z_MIN
                            new_data[i, z_min_lb_f] = c_data[tag, 16]
                            # Z_MAX
                            new_data[i, z_max_lb_f] = c_data[tag, 17]
                            # ODDS in redshift fitting
                            new_data[i, odds_lb_f] = c_data[tag, 18]
                            # # FITCLASS
                            # new_data[i, fitclass_lb] = c_data[tag, 7]
                            # # MASK
                            # new_data[i, mask_lb] = c_data[tag, 9]
                            # # weight
                            # new_data[i, weight_lb] = c_data[tag, 6]
                            # # e1
                            # new_data[i, e1_lb] = c_data[tag, 4]
                            # # e2
                            # new_data[i, e2_lb] = c_data[tag, 5]
                            # # m
                            # new_data[i, m_lb] = c_data[tag, 11]
                            # # c2
                            # new_data[i, c_lb] = c_data[tag, 12]

                            mask[i] += 1

                    idx_mask = mask == 1
                    pairs = idx_mask.sum()
                    idx_mask_2 = mask > 1
                    if pairs == 0:
                        print("%s These two catalogs can't match!!"%field_name)
                        exit(0)
                    elif pairs < f_sp[0]:
                        print("%s Some sources are missing!!"%field_name)
                    if idx_mask_2.sum() > 0:
                        print("%s Overlap!!" % field_name)

                    if field_count == 0:
                        cata_data = new_data
                    else:
                        cata_data = numpy.row_stack((cata_data, new_data))
                    field_count += 1
                except:
                    print(rank, "%s.dat doesn't exist"%field_name)

        # in case of that some thread gets nothing from the file (non-existed)
        if cata_data is not None:
            data_sp = cata_data.shape
        else:
            cata_data = numpy.zeros((1,1))
            data_sp = cata_data.shape

        data_sps = comm.gather(data_sp, root=0)
        num.append(data_sp)
        # npz_name = data_path+"rank_%d_%d.npz"%(i, rank)
        # numpy.savez(npz_name, cata_data)

        if rank > 0:
            comm.Send([cata_data, MPI.DOUBLE], dest=0, tag=rank)
        else:
            if data_sp[0] > 1 and data_sp[1] > 1:
                stack_pool = [cata_data]
            else:
                stack_pool = []

            for procs in range(1, cpus):
                recvs = numpy.empty(data_sps[procs], dtype=numpy.double)
                comm.Recv(recvs, source=procs, tag=procs)
                if data_sps[procs][0] > 1 and data_sps[procs][1] > 1:
                    stack_pool.append(recvs)

            for stack_count, arr in enumerate(stack_pool):
                if stack_count == 0:
                    recv_buffer = arr
                else:
                    recv_buffer = numpy.row_stack((recv_buffer, arr))

            h5f = h5py.File(h5f_path,"r+")

            # each shear estimator is 2 times the one in the paper,
            # Zhang et al. 2017 ApJ 834,
            # G1 = 2*G1_paper, G2 = 2*G2_paper
            # And N = 2N_paper, U = -2U_paper, V = -2V_paper
            recv_buffer[:,mg1_lb] = recv_buffer[:,mg1_lb] / 2
            recv_buffer[:,mg2_lb] = recv_buffer[:,mg2_lb] / 2
            recv_buffer[:,mn_lb] = recv_buffer[:,mn_lb] / 2
            recv_buffer[:,mu_lb] = -recv_buffer[:,mu_lb] / 2
            recv_buffer[:,mv_lb] = -recv_buffer[:,mv_lb] / 2

            h5f["/w_%d" % area_id] = recv_buffer
            h5f.close()
            print("RECV: ", recv_buffer.shape)
        comm.Barrier()
    print(rank, "the galaxy number in each area", num)
    # stack the sub-catalogs from each area
    if rank == 0:
        h5f = h5py.File(h5f_path,"r+")
        for area_id in range(1,area_num+1):
            temp_s = h5f["/w_%d" % area_id].value
            if area_id == 1:
                data = temp_s
            else:
                data = numpy.row_stack((data, temp_s))
            print("Totally, %d galaxies are detected in W_%d" % (len(temp_s), area_id))

        h5f["/total"] = data
        h5f.close()

    t2 = time.time()
    if rank == 0:
        print("%.2f, %d"%(t2-t1, len(data)))
############################# Fourier_quad data collection #######################################


############################# CFHTLenS catalog collection ########################################
# combine the data of each field into one big catalog
if cmd == "collect_cfht":

    t1 = time.time()

    logger.info("RANK: %03d. Begin..."%rank)

    h5f_path = data_path + "cfht_cata.hdf5"
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

            if os.path.exists(c_data_path):
                # try:

                # c_data = numpy.loadtxt(c_data_path)
                h5f_c = h5py.File(c_data_path,"r")
                c_data = h5f_c["/data"].value
                h5f_c.close()

                if field_count == 0:
                    cata_data = c_data

                else:
                    cata_data = numpy.row_stack((cata_data, c_data))
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


############################# Fourier_quad data cutoff ###########################################
if cmd == "select_fourier":
    t1 = time.time()

    h5f_path = data_path + "fourier_cata_%s.hdf5" % result_source
    h5f_path_cut = data_path + "fourier_cata_%s_cut.hdf5" % result_source

    if rank == 0:
        h5f = h5py.File(h5f_path_cut, "w")
        h5f.close()
    comm.Barrier()

    if rank < area_num:
        h5f = h5py.File(h5f_path, "r")
        cata_data = h5f["/w_%d"%(rank+1)].value
        h5f.close()

        # cut off
        flux_alt_idx = cata_data[:, flux_alt_lb] >= flux_alt_thresh
        nstar_idx = cata_data[:, nstar_lb] >= nstar_thresh
        total_area_idx = cata_data[:, total_area_lb] >= total_area_thresh

        fg1 = numpy.abs(cata_data[:, field_g1_lb])
        fg2 = numpy.abs(cata_data[:, field_g2_lb])

        fg1_idx = fg1 <= field_g1_bound
        fg2_idx = fg2 <= field_g2_bound

        idxz_1 = cata_data[:, z_lb_f] <= z_max
        idxz_2 = cata_data[:, z_lb_f] >= z_min

        cut_idx = flux_alt_idx & nstar_idx & total_area_idx & fg1_idx & fg2_idx & idxz_1 & idxz_2

        fg1 = fg1[cut_idx]
        fg2 = fg2[cut_idx]

        # the esitmators
        mn = cata_data[:, mn_lb][cut_idx]
        mu = cata_data[:, mu_lb][cut_idx]
        mv = cata_data[:, mv_lb][cut_idx]
        # correction due to additive bias and field distortion
        mg1 = cata_data[:, mg1_lb][cut_idx] - (fg1 + c1_correction) * (mn + mu) - fg2 * mv
        mg2 = cata_data[:, mg2_lb][cut_idx] - (fg2 + c2_correction) * (mn - mu) - fg1 * mv

        ra = cata_data[:, ra_lb][cut_idx]
        dec = cata_data[:, dec_lb][cut_idx]
        cos_dec = numpy.abs(numpy.cos((dec/180*numpy.pi)))

        redshift = cata_data[:, z_lb_f][cut_idx]
        z_min = cata_data[:,z_min_lb_f][cut_idx]
        z_max = cata_data[:,z_max_lb_f][cut_idx]
        odds = cata_data[:,odds_lb_f][cut_idx]

        mag = cata_data[:, mag_lb_f][cut_idx]
        starflag = cata_data[:, starflag_lb_f][cut_idx]

        # e1 = cata_data[:,e1_lb][cut_idx]
        # e2 = cata_data[:,e2_lb][cut_idx]
        # weight = cata_data[:,weight_lb][cut_idx]
        # fitclass = cata_data[:,fitclass_lb][cut_idx]
        # mask = cata_data[:,mask_lb][cut_idx]
        # m_bias = cata_data[:,m_lb][cut_idx]
        # c_bias = cata_data[:,c_lb][cut_idx]

        names = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V", "MAG", "COS_DEC",
                 "STARGLAG", "Z_MIN", "Z_MAX", "ODDS"]
            # , "E1", "E2", "WEIGHT", "FITCLASS",  "MASK", "M", "C", ]
        datas = [redshift, ra, dec, mg1, mg2, mn, mu, mv, mag, cos_dec, starflag, z_min, z_max, odds]
                 # e1, e2, weight, fitclass, mask, m_bias, c_bias]
        data_num = len(redshift)

    comm.Barrier()

    for area_id in range(area_num):
        if rank == area_id:
            h5f = h5py.File(h5f_path_cut, "r+")
            for i in range(len(names)):
                h5f["/w_%d/%s"%(rank+1,names[i])] = datas[i]
            h5f.close()
        comm.Barrier()
    t2 = time.time()
    if rank == 0:
        print(t2-t1, data_num)
############################# Fourier_quad data cutoff ###########################################


############################# CFHTLenS catlog cutoff #############################################
if cmd == "select_cfht":
    t1 = time.time()

    h5f_path = data_path + "cfht_cata.hdf5"
    h5f_path_cut = data_path + "cfht_cata_cut.hdf5"
    if rank == 0:
        h5f = h5py.File(h5f_path_cut, "w")
        h5f.close()
    comm.Barrier()

    if rank < area_num:
        h5f = h5py.File(h5f_path, "r")
        cata_data = h5f["/w_%d"%(rank+1)].value
        h5f.close()

        # cut off
        epsilon = 0.00001
        # weight > 0
        idx_weight = cata_data[:, weight_lb_c] > 0 - epsilon
        # mask <= 1
        idx_mask = cata_data[:, mask_lb_c] <= 1 + epsilon
        # FITCLASS = 0
        idx_fitclass = numpy.abs(cata_data[:, fitclass_lb_c]) < 0 + epsilon
        # MAG < 24.7
        idx_mag = cata_data[:, mag_lb_c] < 24.7
        # Redshift
        idxz_1 = cata_data[:, z_lb_c] <= z_max
        idxz_2 = cata_data[:, z_lb_c] >= z_min
        idxz_3 = cata_data[:, odds_lb_c] > 0.5 - epsilon

        cut_idx = idx_weight & idx_mask & idx_fitclass & idx_mag & idxz_1 & idxz_2 & idxz_3

        ra = cata_data[:, ra_lb_c][cut_idx]
        dec = cata_data[:, dec_lb_c][cut_idx]
        cos_dec = numpy.abs(numpy.cos((dec/180*numpy.pi)))

        redshift = cata_data[:, z_lb_c][cut_idx]
        z_min = cata_data[:,z_min_lb_c][cut_idx]
        z_max = cata_data[:,z_max_lb_c][cut_idx]
        odds = cata_data[:,odds_lb_c][cut_idx]

        mag = cata_data[:, mag_lb_c][cut_idx]
        starflag = cata_data[:, starflag_lb_c][cut_idx]

        e1 = cata_data[:,e1_lb_c][cut_idx]
        e2 = cata_data[:,e2_lb_c][cut_idx]
        weight = cata_data[:,weight_lb_c][cut_idx]

        fitclass = cata_data[:,fitclass_lb_c][cut_idx]

        mask = cata_data[:,mask_lb_c][cut_idx]

        m_bias = cata_data[:,m_lb_c][cut_idx]
        c_bias = cata_data[:,c_lb_c][cut_idx]

        names = ["Z", "RA", "DEC", "MAG", "COS_DEC", "STARGLAG", "Z_MIN", "Z_MAX", "ODDS",
                 "E1", "E2", "WEIGHT", "FITCLASS", "MASK", "M", "C"]
        datas = [redshift, ra, dec, mag, cos_dec, starflag, z_min, z_max, odds,
                 e1, e2, weight, fitclass, mask, m_bias, c_bias]
        data_num = len(redshift)

    comm.Barrier()

    for area_id in range(area_num):
        if rank == area_id:
            h5f = h5py.File(h5f_path_cut, "r+")
            for i in range(len(names)):
                h5f["/w_%d/%s"%(rank+1,names[i])] = datas[i]
            h5f.close()
        comm.Barrier()
    t2 = time.time()
    if rank == 0:
        print(t2-t1, data_num)
############################# CFHTLenS catlog cutoff #############################################