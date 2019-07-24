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


# The new Fourier_Quad catalog differs from the old version!!!
# collect: collect the data from the files of each field. It creates the "fourier_cata.hdf5" in
#           the parent directory of the one contain the field catalog.
#           If the catalog file doesn't exist, run it firstly !!!.
#           It will add the redshift parameters from CFHT catalog into the finial catalog.

# select: select the galaxy

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

envs_path = "%s/work/envs/envs.dat"%my_home

gets_item = [["cfht", "cfht_path_catalog", "0"], ["gg_lensing", "ggl_path_data", "0"]]
path_items = tool_box.config(envs_path, ["get", "get"], gets_item)

cata_path, data_path = path_items

cfht_cata_path = cata_path + "cfht_cata/"
fourier_cata_path = cata_path + "fourier_cata_new/"


gets_item = [["fresh_para_idx", "nstar", "0"], ["fresh_para_idx", "flux_alt", "0"],
             ["fresh_para_idx", "ra", "0"], ["fresh_para_idx", "dec", "0"],
             ["fresh_para_idx", "gf1", "0"], ["fresh_para_idx", "gf2", "0"],
             ["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
             ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
             ["fresh_para_idx", "h2", "0"], ["fresh_para_idx", "total_area", "0"]]
gets = ["get" for i in range(len(gets_item))]
para_items = tool_box.config(envs_path, gets, gets_item)

############################# Fourier_Quad Option ################################################
# the column label in Fourier_Quad (new) catalog,
# used in Fourier_quad related option


file_header = "RA  DEC ig  g1  g2  de  h1  h2  Z   Mag_i   FLUX_F  nstar   Z_MIN   ZMAX    ODDS    PZ1  PZ2"
# the label in Fourier_Quad catalog
ra_lb = 0
dec_lb = 1

ig_lb = 2

mg1_lb = 3
mg2_lb = 4
mn_lb = 5
mu_lb = 6
mv_lb = 7

z_lb = 8
mag_lb = 9
flux_alt_lb = 10
nstar_lb = 11

# it would change as more parameters added into the file
data_col = 12
add_col = 5
# the added columns
z_min_lb = 12
z_max_lb = 13
odds_lb = 14
pz1_lb = 15
pz2_lb = 16

# the labels of the added parameters in the CFHT catalog
ra_lb_c = 0
dec_lb_c = 1
z_lb_c = 10
z_min_lb_c = 16
z_max_lb_c = 17
odds_lb_c = 18


# cut off for the "select" option
flux_alt_thresh = 4
nstar_thresh = 12
z_min, z_max = 0.0, 15
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
            f_data_path = fourier_cata_path + "%s/%s_shear.cat"%(field_name, field_name)
            f_data_path_new = fourier_cata_path + "%s/%s_shear_new.cat"%(field_name, field_name)
            # c_data_path = cfht_cata_path + "%s_new.dat"%field_name
            c_data_path = cfht_cata_path + "field_dat/%s.hdf5"%field_name
            pz_data_path = cfht_cata_path + "field_dat/%s-pz.dat"%field_name

            print(f_data_path, c_data_path)

            if os.path.exists(f_data_path) and os.path.exists(c_data_path):
                # try:
                logger.info("Begin to match...")
                # the file in Fourier_Quad catalog may be empty
                f_data = numpy.loadtxt(f_data_path)
                f_sp = f_data.shape
                # the additional parameters, if non-exist in the cfht catalog, labeled by -1
                new_data = numpy.zeros((f_sp[0], data_col + add_col), dtype=numpy.double) - 1
                new_data[:,:data_col] = f_data[:,:data_col]

                # c_data = numpy.loadtxt(c_data_path)
                pz_data = numpy.loadtxt(pz_data_path)

                h5f_c = h5py.File(c_data_path,"r")
                c_data = h5f_c["/data"].value
                h5f_c.close()

                c_sp = c_data.shape
                print(c_sp)
                check_buf = numpy.zeros((f_sp[0], 2)) - 99

                # f_sp[0] > c_sp[0]
                pairs = 0
                mask = numpy.zeros((f_sp[0],), dtype=numpy.intc)
                for i in range(f_sp[0]):
                    # fortran start from 1 not 0
                    # the row labels corresponding to the galaxies in cfht catalog
                    # every source in Fourier_Quad ext_cata has label corresponding
                    # to the row in the CFHTLenS catalog
                    tag = int(f_data[i, ig_lb] - 1)
                    # check
                    del_radius = numpy.abs(f_data[i, ra_lb] - c_data[tag, ra_lb_c]) + \
                                 numpy.abs(f_data[i, dec_lb] - c_data[tag, dec_lb_c])
                    del_z = numpy.abs(f_data[i, z_lb] - c_data[tag, z_lb_c])

                    check_buf[i, 0] = del_radius
                    check_buf[i, 1] = del_z

                    if del_radius < 0.00005 and del_z < 0.00005:
                        # Z_MIN
                        new_data[i, z_min_lb] = c_data[tag, z_min_lb_c]
                        # Z_MAX
                        new_data[i, z_max_lb] = c_data[tag, z_max_lb_c]
                        # ODDS in redshift fitting
                        new_data[i, odds_lb] = c_data[tag, odds_lb_c]
                        # PZ1 from Dong FY about the Z quality
                        new_data[i, pz1_lb] = pz_data[tag, 0]
                        # PZ2 from Dong FY about the Z quality
                        new_data[i, pz2_lb] = pz_data[tag, 1]

                        mask[i] += 1
                numpy.savetxt(fname=f_data_path_new, X=new_data, header=file_header)
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
                logger.info("write file %s. MASK: max = %d, min = %d. Check buff: max = %.6f" % (f_data_path_new, mask.max(), mask.min(), check_buf.max()))
                # except:
                #     print(rank, "%s.cat doesn't exist"%field_name)

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

    h5f_path = fourier_cata_path + "fourier_cata.hdf5"
    h5f_path_cut = data_path + "fourier_cata_cut.hdf5"

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

        idxz_1 = cata_data[:, z_lb] <= z_max
        idxz_2 = cata_data[:, z_lb] >= z_min

        cut_idx = flux_alt_idx & nstar_idx & idxz_1 & idxz_2

        # the esitmators
        mn = cata_data[:, mn_lb][cut_idx]
        mu = cata_data[:, mu_lb][cut_idx]
        mv = cata_data[:, mv_lb][cut_idx]
        # correction due to additive bias and field distortion
        mg1 = cata_data[:, mg1_lb][cut_idx]
        mg2 = cata_data[:, mg2_lb][cut_idx]

        ra = cata_data[:, ra_lb][cut_idx]
        dec = cata_data[:, dec_lb][cut_idx]
        cos_dec = numpy.abs(numpy.cos((dec/180*numpy.pi)))

        redshift = cata_data[:, z_lb][cut_idx]
        z_min = cata_data[:,z_min_lb][cut_idx]
        z_max = cata_data[:,z_max_lb][cut_idx]
        odds = cata_data[:,odds_lb][cut_idx]
        pz1 = cata_data[:, pz1_lb][cut_idx]
        pz2 = cata_data[:, pz2_lb][cut_idx]

        mag = cata_data[:, mag_lb][cut_idx]

        names = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V", "MAG", "COS_DEC",
                 "Z_MIN", "Z_MAX", "ODDS", "PZ1", "PZ2"]

        datas = [redshift, ra, dec, mg1, mg2, mn, mu, mv, mag, cos_dec, z_min, z_max, odds, pz1, pz2]

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


