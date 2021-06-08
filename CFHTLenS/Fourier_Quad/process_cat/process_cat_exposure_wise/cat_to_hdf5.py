import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
import warnings

warnings.filterwarnings('error')


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

total_path = argv[1]
mode = argv[2]

if mode == "cata_name":
    if rank == 0:

        # with open("anomaly_expo.dat", "r") as f:
        #     contents = f.readlines()
        # anomaly_expo = []
        # anomaly_val = []
        # for cc in contents:
        #     val, field_nm, expo = cc.rstrip().split()[1:4]
        #     anomaly_expo.append(int(expo.split("p")[0]))
        #     anomaly_val.append(float(val))
        #     # print(anomaly_expo[-1],anomaly_val[-1])

        files_nm = os.listdir(total_path + "/cat_ori")
        all_files = []

        for fnm in files_nm:
            if "_all.cat" in fnm:
                all_files.append("%s\n"%fnm.split(".")[0])

        with open(total_path + "/cat_inform/exposure_name.dat", "w") as f:
            f.writelines(all_files)

        print(len(all_files), " exposures")

elif mode == "hdf5_cata":
    # convert the .dat to .hdf5
    # CFHT_pz_path = argv[3]
    # # read the photoz data
    # # RA DEC Z_B Z_B_MIN Z_B_MAX ODDS Z_E Pz_full(70 cols)
    # for i in range(numprocs):
    #     if i == rank:
    #         h5f = h5py.File(CFHT_pz_path,"r")
    #         pz_data = h5f["/data"][()]
    #         h5f.close()
    #         # print("%d Read Pz cata"%rank)
    #
    #     comm.Barrier()
    # comm.Barrier()

    # zbin = numpy.array([0.025 + i * 0.05 for i in range(70)])

    exposures_candidates = []
    with open(total_path + "/cat_inform/exposure_name.dat", "r") as f:
        contents = f.readlines()
    for expo_name in contents:
        exposures_candidates.append(expo_name.split("\n")[0])

    exposures_candidates_sub = tool_box.alloc(exposures_candidates, numprocs, "seq")[rank]

    exposures_candidates_avail_sub = []
    exception_sub = []

    for fns in exposures_candidates_sub:
        # read the the field data
        expo_src_path = total_path + "/cat_ori/%s.cat" % fns
        expo_h5_path = total_path + "/cat_hdf5/%s.hdf5" % fns
        try:
            src_data = numpy.loadtxt(expo_src_path, dtype=numpy.float32)
            row, col = src_data.shape
        except:
            if os.path.exists(expo_src_path):
                log_inform = "%d Failed in reading %s %d Bytes !\n" % (
                    rank, expo_src_path, os.path.getsize(expo_src_path))
            else:
                log_inform = "%d can't find %s!\n" % (rank, expo_src_path)
            exception_sub.append(log_inform)
            row, col = 0, 0

        if row > 0:
            # # for Z_B_MIN, Z_B_MAX, ODDS, Z_E, sum(Pz)
            # dst_data = numpy.zeros((row, col+5), dtype=numpy.float32)
            # dst_data[:,:col] = src_data
            #
            # ra_min, ra_max = src_data[:,0].min(), src_data[:,0].max()
            # dra = (ra_max - ra_min)*0.06
            # dec_min, dec_max = src_data[:,1].min(), src_data[:,1].max()
            # ddec = (dec_max - dec_min) * 0.06
            #
            # idx1 = pz_data[:,0] >= ra_min - dra
            # idx2 = pz_data[:,0] <= ra_max + dra
            # idx3 = pz_data[:,1] >= dec_min - ddec
            # idx4 = pz_data[:,1] <= dec_max + ddec
            # idx = idx1 & idx2 & idx3 & idx4
            # pz_data_sub = pz_data[idx]
            # pz_data_sub_num = idx.sum()
            # labels = numpy.arange(0, pz_data_sub_num)
            #
            # src_num = src_data.shape[0]
            # for i in range(src_num):
            #     diff_rad = numpy.abs(pz_data_sub[:,0] - src_data[i,0]) + numpy.abs(pz_data_sub[:,1] - src_data[i,1])
            #     diff_z = numpy.abs(pz_data_sub[:,2] - src_data[i,10])
            #
            #     diff_rad_min = diff_rad.min()
            #     diff_z_min = diff_z.min()
            #
            #     idx_1 = diff_rad == diff_rad.min()
            #     idx_2 = diff_z == diff_z.min()
            #     idx_ = idx_1 & idx_2
            #     match_num = idx_.sum()
            #
            #     if diff_rad_min >= 0.0001:
            #         print("Can't find match points!!! diff_radius_min: %.5f"%diff_rad_min)
            #         print(labels[idx_1])
            #         exit()
            #
            #     if match_num == 1:
            #         target_idx = labels[idx_]
            #     else:
            #         print("Find %d pts!!!"%match_num)
            #         print(idx_1.sum(), idx_2.sum(), idx_.sum())
            #         print(diff_rad.min(),diff_z.min())
            #         print("%.4f %.4f  %.3f"%(src_data[i,0],src_data[i,1], src_data[i,10]))
            #         exit()
            #
            #     # Z_B_MIN
            #     dst_data[i,col] = pz_data_sub[target_idx, 3]
            #     # Z_B_MAX
            #     dst_data[i,col+1] = pz_data_sub[target_idx, 4]
            #     # ODDS
            #     dst_data[i,col+2] = pz_data_sub[target_idx, 5]
            #     # expectation of Z
            #     dst_data[i,col+3] = pz_data_sub[target_idx, 6]
            #     # sum(Pz)
            #     dst_data[i,col+4] = pz_data_sub[target_idx, 7]

            # Nan check
            idx = numpy.isnan(src_data)
            if idx.sum() > 0:
                print("Find Nan in ", expo_src_path)

            h5f_expo = h5py.File(expo_h5_path, "w")
            h5f_expo["/data"] = src_data
            h5f_expo.close()

            exposures_candidates_avail_sub.append(expo_h5_path + "\n")

    exposures_candidates_avail = comm.gather(exposures_candidates_avail_sub, root=0)
    exception_collection = comm.gather(exception_sub, root=0)

    comm.Barrier()
    if rank == 0:
        exception_all = []
        for ec in exception_collection:
            exception_all.extend(ec)

        exposures_avail = []
        for fsb in exposures_candidates_avail:
            exposures_avail.extend(fsb)
        with open(total_path + "/cat_inform/exposure_avail.dat","w") as f:
            f.writelines(exposures_avail)

        exception_all.append("Totally: %d/%d available exposures\n"
                             "The expectation of Z calculated from the P(z) and normalized P(z), "
                             "and the ODDS have been added to the hdf5 cat\n"%(len(exposures_avail),len(exposures_candidates)))
        with open("log.dat","w") as f:
            f.writelines(exception_all)


elif mode == "field_expo":

    if rank == 0:
        field_expos = []
        with open(total_path + "/cat_inform/field_expo.dat", "r") as f:
            cc = f.readlines()
        expo_num = 0
        for c in cc:
            nms = c.split("\n")[0].split("\t")
            temp = [nms[0]]
            for expo_name in nms[1:]:
                expo_path = total_path + "/cat_hdf5/%s_all.hdf5"%expo_name
                if os.path.exists(expo_path):
                    temp.append(expo_name)
                    expo_num += 1
                else:
                    print("Not found %s"%expo_path)
            temp_str = "\t".join(temp) + "\n"
            field_expos.append(temp_str)
        print(expo_num)
        with open(total_path + "/cat_inform/field_expo_avail.dat", "w") as f:
            f.writelines(field_expos)
    comm.Barrier()


elif mode == "position":
    expo_files = []
    with open(total_path + "/cat_inform/field_expo_avail.dat", "r") as f:
        cc = f.readlines()
    field_num = len(cc)

    max_expo_num = 0
    for c in cc:
        nms = c.split("\n")[0].split("\t")
        expo_num = len(nms) - 1
        if expo_num > max_expo_num:
            max_expo_num = expo_num
        expo_files.append(nms[1:])

    expo_pos = numpy.zeros((field_num, max_expo_num, 8), dtype=numpy.float32)

    field_label = [i for i in range(field_num)]
    field_sub = tool_box.alloc(field_label, numprocs,"seq")[rank]

    for i in field_sub:
        for j, expo_name in enumerate(expo_files[i]):
            file_path = total_path + "/cat_hdf5/%s_all.hdf5"%expo_name

            h5f = h5py.File(file_path, "r")
            data = h5f["/data"][()]
            h5f.close()

            ra, dec = data[:,0], data[:,1]
            ra_min, ra_max = ra.min(), ra.max()
            dec_min, dec_max = dec.min(), dec.max()

            dra = (ra_max - ra_min)/2
            ddec = (dec_max - dec_min)/2

            ra_cent = (ra_max + ra_min)/2
            dec_cent = (dec_max + dec_min)/2

            expo_pos[i][j] = ra_cent, dec_cent, ra_min, ra_max, dra, dec_min, dec_max, ddec

    comm.Barrier()

    if rank > 0:
        comm.Send([expo_pos, MPI.FLOAT], dest=0, tag=rank)
    else:
        for procs in range(1, numprocs):
            recvs = numpy.empty((field_num, max_expo_num, 8), dtype=numpy.float32)
            comm.Recv(recvs, source=procs, tag=procs)
            expo_pos += recvs

        h5f = h5py.File(total_path + "/cat_inform/expo_pos.hdf5", "w")
        h5f["/data"] = expo_pos
        h5f.close()

    comm.Barrier()


else:
    # collection

    source_list_nm = argv[2]
    result_nm = argv[3]

    expos = []
    gal_num = []
    with open(total_path + "/"+source_list_nm, "r") as f:
        conts = f.readlines()
    for nm in conts:
        informs = nm.split("\n")
        expos.append(informs[0])
        gal_num.append(informs[2])

    if rank == 0:
        print(len(expos)," exposures")

    my_sub_area_list = tool_box.alloc(expos,numprocs)[rank]
    # print(rank, i, len(my_sub_area_list))

    if len(my_sub_area_list) > 0:
        for tag, expo_path in enumerate(my_sub_area_list):

            h5f = h5py.File(expo_path,"r")
            temp = h5f["/data"][()]

            # Nan check
            idx = numpy.isnan(temp)
            if idx.sum() > 0:
                print("Find Nan ", expo_path)
            if tag == 0:
                data = temp
            else:
                data = numpy.row_stack((data, temp))
            h5f.close()

        sp = data.shape
    else:
        sp = (0,0)

    sp_total = comm.gather(sp, root=0)
    # print(i,rank, data.shape, data.dtype, sp, data[0,:5])
    comm.Barrier()

    if rank > 0 and sp[0] > 0:
        comm.Send([data,MPI.FLOAT], dest=0, tag=rank)
    else:
        for ir in range(1, numprocs):
            if sp_total[ir][0] > 0:
                recv_buf = numpy.empty(sp_total[ir],dtype=numpy.float32)
                comm.Recv(recv_buf,source=ir, tag=ir)
                data = numpy.row_stack((data, recv_buf))

        # Nan check
        idx = numpy.isnan(data)
        if idx.sum() > 0:
            print("Find Nan in final data")

        idx_v = numpy.invert(idx)
        h5f = h5py.File(total_path + "/%s"%result_nm, "w")
        h5f["/data"] = data
        h5f.close()
        print("Totally %d (%d) galaxies"%(data.shape[0], sum(gal_num)))


    comm.Barrier()