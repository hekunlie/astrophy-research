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
import matplotlib.pyplot as plt
import numpy
import time


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]
cmds = ["collect", "select"]

# collect: collect the data from the files. run it firstly.
#           add MAG, Z, Z_MIN, Z_MAX, ODDS, starflag, weight , and mask
#           to end of each row in Fourier catalog
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
odds_lb = -1
z_max_lb = -2
z_min_lb = -3
z_lb = -4
mag_lb = -5
starflag_lb = -6

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


# data collection
if cmd == "collect":

    t1 = time.time()

    h5f_path = data_path + "cata_%s.hdf5" % result_source
    if rank == 0:
        h5f = h5py.File(h5f_path, "w")
        h5f.close()
        buffer = []

    c_dicts, c_fields = tool_box.field_dict(cfht_cata_path+"nname.dat")
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
                    new_data = numpy.zeros((f_sp[0], f_sp[1] + 6), dtype=numpy.double) - 1
                    new_data[:,:f_sp[1]] = f_data

                    # c_data = numpy.loadtxt(c_data_path)
                    h5f_c = h5py.File(c_data_path,"r")
                    c_data = h5f_c["/data"].value
                    h5f_c.close()

                    c_sp = c_data.shape
                    # f_sp[0] > c_sp[0]
                    pairs = 0
                    for i in range(f_sp[0]):
                        # fortran start from 1 not 0
                        # the row labels corresponding to the galaxies in cfht catalog
                        # every source in Fourier_Quad ext_cata has label corresponding
                        # to the row in the CFHTLenS catalog
                        tag = int(f_data[i, 0] - 1)
                        # check
                        del_radius = numpy.abs(f_data[i,ra_lb] - c_data[tag, 0]) + numpy.abs(f_data[i,dec_lb] - c_data[tag, 1])

                        if del_radius <= 0.000005:
                            # starflag
                            new_data[i, f_sp[1]] = c_data[tag, 14]
                            # magnitude
                            new_data[i, f_sp[1] + 1] = c_data[tag, 15]
                            # redshift
                            new_data[i, f_sp[1] + 2] = c_data[tag, 10]
                            # Z_MIN
                            new_data[i, f_sp[1] + 3] = c_data[tag, 16]
                            # Z_MAX
                            new_data[i, f_sp[1] + 4] = c_data[tag, 17]
                            # ODDS in redshift fitting
                            new_data[i, f_sp[1] + 5] = c_data[tag, 18]
                            pairs += 1

                    if pairs == 0:
                        print("These two catalogs can't match!!")
                        exit(0)

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
            cata_data = numpy.array([[1]])
            data_sp = cata_data.shape

        data_sps = comm.gather(data_sp, root=0)
        num.append(data_sp)
        # npz_name = data_path+"rank_%d_%d.npz"%(i, rank)
        # numpy.savez(npz_name, cata_data)

        if rank > 0:
            comm.Send(cata_data, dest=0, tag=rank)
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

            h5f = h5py.File(h5f_path)

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
        # fig1 = plt.figure(figsize=(14, 14))
        # fig2 = plt.figure(figsize=(14, 14))
        # fig3 = plt.figure(figsize=(14, 14))

        h5f = h5py.File(h5f_path,"r+")

        for area_id in range(1,area_num+1):
            temp_s = h5f["/w_%d" % area_id].value

            # # select the data as CFTHLenS
            # red_z = temp[:, -3]
            # idxz_1 = red_z < z_max
            # idxz_2 = red_z > z_min
            # temp_s = temp[idxz_1&idxz_2]
            #
            # # show the range of RA & DEC
            # print(area_id, temp_s.shape,temp_s[:,12].min()*60, temp_s[:,12].max()*60,
            #       temp_s[:,13].min()*60, temp_s[:,13].max()*60)
            # ax1 = fig1.add_subplot(220 + area_id)
            # ax1.scatter(temp_s[:,12]*60, temp_s[:,13]*60, s=0.3)
            #
            # # histogram of redshift
            # ax2 = fig2.add_subplot(220 + area_id)
            # ax2.hist(red_z, 50)
            #
            # # of magnitude
            # mag = temp_s[:, -2]
            # idxm_1 = mag > 0
            # idxm_2 = mag < 30
            # ax3 = fig3.add_subplot(220 + area_id)
            # ax3.hist(mag[idxm_1&idxm_2], 50)

            if area_id == 1:
                data = temp_s
            else:
                data = numpy.row_stack((data, temp_s))
            print("Totally, %d galaxies are detected in W_%d" % (len(temp_s), area_id))

        # fig1.savefig(data_path + "pic/Ra_dec_%s.png"%result_source)
        # fig2.savefig(data_path + "pic/Z_%s.png"%result_source)
        # fig3.savefig(data_path + "pic/Mag_%s.png"%result_source)
        plt.close()

        h5f["/total"] = data
        h5f.close()

    t2 = time.time()
    if rank == 0:
        print(t2-t1, len(data))

if cmd == "select":
    t1 = time.time()

    h5f_path = data_path + "cata_%s.hdf5" % result_source
    h5f_path_cut = data_path + "cata_%s_cut.hdf5" % result_source
    h5f_red_dist_path = data_path + "redshift.hdf5"
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

        idxz_1 = cata_data[:, z_lb] <= z_max
        idxz_2 = cata_data[:, z_lb] >= z_min

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

        redshift = cata_data[:, z_lb][cut_idx]
        z_min = cata_data[:,z_min_lb][cut_idx]
        z_max = cata_data[:,z_max_lb][cut_idx]
        odds = cata_data[:,odds_lb][cut_idx]

        mag = cata_data[:, mag_lb][cut_idx]
        starflag = cata_data[:, starflag_lb][cut_idx]

        gal_num = len(mg1)
        ra_min, ra_max = ra.min(),ra.max()
        dec_min, dec_max = dec.min(), dec.max()

        # distance = numpy.zeros_like(redshift)
        # for i in range(gal_num):
        #     tag = tool_box.find_near(redshift_refer, redshift[i])
        #     distance[i] = distance_refer[tag]
        names = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V", "MAG", "COS_DEC", "STARGLAG","Z_MIN", "Z_MAX","ODDS"]
        datas = [redshift, ra, dec, mg1, mg2, mn, mu, mv, mag, cos_dec, starflag, z_min, z_max, odds]
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
