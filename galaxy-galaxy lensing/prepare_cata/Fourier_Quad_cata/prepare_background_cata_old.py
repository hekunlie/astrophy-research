import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
from Fourier_Quad import Fourier_Quad
import tool_box
import h5py
from mpi4py import MPI
from sys import argv
import numpy
import time
from subprocess import Popen


# The new Fourier_Quad catalog differs from the old version!!!
# collect: collect the data from the files of each field. It creates the "fourier_cata.hdf5" in
#           the parent directory of the one contain the field catalog.
#           If the catalog file doesn't exist, run it firstly !!!.
#           It will add the redshift parameters from CFHT catalog into the finial catalog.

# select: select the galaxy
# grid:

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]
cmds = ["collect", "select", "grid"]

if cmd not in cmds:
    if rank == 0:
        print("parameter must be one of ", cmds)
    exit()

area_num = 4

cata_path = "/mnt/perc/hklee/CFHT/catalog/"
data_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/"

cfht_cata_path = cata_path + "cfht_cata/"
fourier_cata_path = cata_path + "fourier_cata_old/"

envs_path = "%s/work/envs/envs.dat"%my_home
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

add_col = 5

mag_lb_f = 21
z_lb_f = 22
z_min_lb_f = 23
z_max_lb_f = 24
odds_lb_f = 25

# the labels of the added parameters in the CFHT catalog
ra_lb_c = 0
dec_lb_c = 1
z_lb_c = 10
mag_lb_c = 15
z_min_lb_c = 16
z_max_lb_c = 17
odds_lb_c = 18

# cut off
flux_alt_thresh = 5.5
nstar_thresh = 12
total_area_thresh = 1
field_g1_bound = 0.005
field_g2_bound = 0.0075
z_min, z_max = 0.0, 15
c1_correction = -0.000578
c2_correction = 0.000493

############################# Fourier_Quad Option ################################################


logger = tool_box.get_logger("%s/work/test/log/%d.dat"%(my_home,rank))

############################# Fourier_quad data collection #######################################
# combine the data of each field into one big catalog
if cmd == "collect":

    t1 = time.time()

    h5f_path = fourier_cata_path + "fourier_cata.hdf5"
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
            f_data_path = fourier_cata_path + "%s/result_ext/%s_shear.dat"%(field_name, field_name)
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
                    check_buf = numpy.zeros((f_sp[0],1)) - 99
                    for i in range(f_sp[0]):
                        # fortran start from 1 not 0
                        # the row labels corresponding to the galaxies in cfht catalog
                        # every source in Fourier_Quad ext_cata has label corresponding
                        # to the row in the CFHTLenS catalog
                        tag = int(f_data[i, 0] - 1)
                        # check
                        del_radius = numpy.abs(f_data[i,ra_lb] - c_data[tag, 0]) + numpy.abs(f_data[i,dec_lb] - c_data[tag, 1])

                        check_buf[i, 0] = del_radius

                        if del_radius < 0.00005:
                            # magnitude
                            new_data[i, mag_lb_f] = c_data[tag, mag_lb_c]
                            # redshift
                            new_data[i, z_lb_f] = c_data[tag, z_lb_c]
                            # Z_MIN
                            new_data[i, z_min_lb_f] = c_data[tag, z_min_lb_c]
                            # Z_MAX
                            new_data[i, z_max_lb_f] = c_data[tag, z_max_lb_c]
                            # ODDS in redshift fitting
                            new_data[i, odds_lb_f] = c_data[tag, odds_lb_c]

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
                    logger.info("MASK: max = %d, min = %d. Check buff: max = %.6f" % (mask.max(), mask.min(),check_buf.max()))
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


############################# Fourier_quad data cutoff ###########################################
if cmd == "select":
    t1 = time.time()

    block_scale = float(argv[2])
    margin = 0.1*block_scale

    mg_bin_num = 8

    h5f_path = fourier_cata_path + "fourier_cata.hdf5"
    h5f_path_cut = data_path + "fourier_cata_old/fourier_cata_cut.hdf5"

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
        # the minus arises from the true RA-axis is in opposite direction
        mg2 = -(cata_data[:, mg2_lb][cut_idx] - (fg2 + c2_correction) * (mn - mu) - fg1 * mv)

        ra = cata_data[:, ra_lb][cut_idx]
        dec = cata_data[:, dec_lb][cut_idx]
        cos_dec = numpy.abs(numpy.cos((dec/180*numpy.pi)))

        redshift = cata_data[:, z_lb_f][cut_idx]
        z_min = cata_data[:,z_min_lb_f][cut_idx]
        z_max = cata_data[:,z_max_lb_f][cut_idx]
        odds = cata_data[:,odds_lb_f][cut_idx]

        mag = cata_data[:, mag_lb_f][cut_idx]

        names = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V", "MAG", "COS_DEC", "Z_MIN", "Z_MAX", "ODDS"]

        datas = [redshift, ra, dec, mg1, mg2, mn, mu, mv, mag, cos_dec, z_min, z_max, odds]

        fq = Fourier_Quad(12,123)
        mg1_bin = fq.set_bin(mg1, mg_bin_num, 10000)
        mg2_bin = fq.set_bin(mg2, mg_bin_num, 10000)

        data_num = len(redshift)

        total_num = ra.shape[0]
        gal_label = numpy.arange(0,data_num)

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

        # galaxy count in each block
        num_in_block = numpy.zeros((1, grid_num), dtype=numpy.intc)
        block_start = numpy.zeros((1, grid_num), dtype=numpy.intc)
        # the galaxy labels in the block
        gal_sequence = numpy.zeros((1, data_num), dtype=numpy.intc)

        for i in range(ny):
            idx_1 = dec >= dec_bin[i]
            idx_2 = dec < dec_bin[i+1]
            idx_ = idx_1 & idx_2
            sub_ra = ra[idx_]
            sub_gal_label = gal_label[idx_]
            ix = i*nx
            for j in range(nx):
                tag = ix + j

                boundy[tag,0] = dec_bin[i]
                boundy[tag,1] = dec_bin[i]
                boundy[tag,2] = dec_bin[i+1]
                boundy[tag,3] = dec_bin[i+1]

                boundx[tag,0] = ra_bin[j]
                boundx[tag,1] = ra_bin[j+1]
                boundx[tag,2] = ra_bin[j]
                boundx[tag,3] = ra_bin[j+1]

                idx_3 = sub_ra >= ra_bin[j]
                idx_4 = sub_ra < ra_bin[j+1]
                idx = idx_3 & idx_4
                num_in_block[0,tag] = idx.sum()
                if tag > 0:
                    block_start[0,tag] = num_in_block[0,:tag-1].sum()
                gal_sequence[0, block_start[0, tag]: block_start[0, tag]+num_in_block[0, tag]] = sub_gal_label[idx]

        print("W%d: %d ==> %d"%(rank, cata_data.shape[0], data_num))
    comm.Barrier()

    for area_id in range(area_num):
        if rank == area_id:
            h5f = h5py.File(h5f_path_cut, "r+")

            h5f["/w_%d/grid_shape" % (rank + 1)] = grid_shape

            h5f["/w_%d/num_in_block" % (rank + 1)] = num_in_block
            h5f["/w_%d/block_start" % (rank + 1)] = block_start
            h5f["/w_%d/gal_in_block" % (rank + 1)] = gal_sequence

            h5f["/w_%d/block_boundy" % (rank + 1)] = boundy
            h5f["/w_%d/block_boundx" % (rank + 1)] = boundx
            h5f["/w_%d/RA_bin" % (rank + 1)] = ra_bin
            h5f["/w_%d/DEC_bin" % (rank + 1)] = dec_bin

            h5f["/w_%d/mg1_bin" % (rank + 1)] = mg1_bin
            h5f["/w_%d/mg2_bin" % (rank + 1)] = mg2_bin

            for i in range(len(names)):
                h5f["/w_%d/%s"%(rank+1,names[i])] = datas[i]

            h5f.close()
        comm.Barrier()
    t2 = time.time()
    comm.Barrier()
    if rank == 0:
        for i in range(area_num):
            cmd = "../add_com_dist %s /w_%d/"%(h5f_path_cut, i+1)
            a = Popen(cmd, shell=True)
            a.wait()
        print("%.2f sec"%(t2-t1))
############################# Fourier_quad data cutoff ###########################################
