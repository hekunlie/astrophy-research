import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
from Fourier_Quad import Fourier_Quad
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
cmds = ["collect","redshift","grid"]

# collect: collect the data from the files. run it firstly
# redshift: divide the data from "collect" into several redshift bins
#           the data will be cut off (nstar, flux_alt, redshift..).
#           run it before "grid"
# grid: divide the data into the grid

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
z_lb = -3

flux_alt_thresh = 6.38
nstar_thresh = 14
total_area_thresh = 1
field_g1_bound = 0.005
field_g2_bound = 0.0075
z_min, z_max = 0.2, 1.3
c1_correction = -0.000287
c2_correction = 0.000609


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
            f_data_path = fourier_cata_path + "%s/%s/%s_shear.dat"%(field_name, result_source, field_name)
            c_data_path = cfht_cata_path + "%s.dat"%field_name

            if os.path.exists(f_data_path) and os.path.exists(c_data_path):
                try:
                    # the file in Fourier_Quad catalog may be empty
                    f_data = numpy.loadtxt(f_data_path)
                    f_sp = f_data.shape
                    # the redshift magnitude and starflag, if non-exist in the cfht catalog, labeled by -1
                    z_mag_col = numpy.zeros((f_sp[0], 3), dtype=numpy.double) - 1
                    f_data = numpy.column_stack((f_data, z_mag_col))

                    c_data = numpy.loadtxt(c_data_path)
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
                        d_ra = abs(f_data[i,ra_lb] - c_data[tag, 0])
                        d_dec = abs(f_data[i,dec_lb] - c_data[tag, 1])

                        if d_ra <= 10**(-4) and d_dec <= 10**(-4):
                            # redshift
                            f_data[i, -3] = c_data[tag, 10]
                            # magnitude
                            f_data[i, -2] = c_data[tag, 15]
                            # starglag
                            f_data[i, -1] = c_data[tag, 14]
                            pairs += 1

                    if pairs == 0:
                        print("These two catalogs can't match!!")
                        exit(0)

                    if field_count == 0:
                        cata_data = f_data
                    else:
                        cata_data = numpy.row_stack((cata_data, f_data))
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
            # And N = -2N_paper, U = -2U_paper, V = -2V_paper
            recv_buffer[:,mg1_lb] = recv_buffer[:,mg1_lb] / 2
            recv_buffer[:,mg2_lb] = recv_buffer[:,mg2_lb] / 2
            recv_buffer[:,mn_lb] = recv_buffer[:,mn_lb] / 2
            recv_buffer[:,mu_lb] = -recv_buffer[:,mu_lb] / 2
            recv_buffer[:,mv_lb] = -recv_buffer[:,mv_lb] / 2

            h5f["/w_%d" % area_id] = recv_buffer
            h5f.close()
            print("RECV: ", recv_buffer.shape)

    print(rank, "the galaxy number in each area", num)
    # stack the sub-catalogs from each area
    if rank == 0:
        fig1 = plt.figure(figsize=(14, 14))
        fig2 = plt.figure(figsize=(14, 14))
        fig3 = plt.figure(figsize=(14, 14))

        h5f = h5py.File(h5f_path,"r+")

        for area_id in range(1,area_num+1):
            temp = h5f["/w_%d" % area_id].value

            # select the data as CFTHLenS
            red_z = temp[:, -3]
            idxz_1 = red_z < z_max
            idxz_2 = red_z > z_min
            temp_s = temp#[idxz_1&idxz_2]

            # show the range of RA & DEC
            print(area_id, temp_s.shape,temp_s[:,12].min()*60, temp_s[:,12].max()*60,
                  temp_s[:,13].min()*60, temp_s[:,13].max()*60)
            ax1 = fig1.add_subplot(220 + area_id)
            ax1.scatter(temp_s[:,12]*60, temp_s[:,13]*60, s=0.3)

            # histogram of redshift
            ax2 = fig2.add_subplot(220 + area_id)
            ax2.hist(red_z, 50)

            # of magnitude
            mag = temp_s[:, -2]
            idxm_1 = mag > 0
            idxm_2 = mag < 30
            ax3 = fig3.add_subplot(220 + area_id)
            ax3.hist(mag[idxm_1&idxm_2], 50)

            if area_id == 1:
                data = temp_s
            else:
                data = numpy.row_stack((data, temp_s))
            print("Totally, %d galaxies are detected in W_%d" % (len(temp_s), area_id))

        fig1.savefig(data_path + "pic/Ra_dec_%s.png"%result_source)
        fig2.savefig(data_path + "pic/Z_%s.png"%result_source)
        fig3.savefig(data_path + "pic/Mag_%s.png"%result_source)
        plt.close()

        h5f["/total"] = data
        h5f.close()

    t2 = time.time()
    if rank == 0:
        print(t2-t1, len(data))