import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
import warnings
from sklearn.cluster import KMeans
import time

warnings.filterwarnings('error')


# parameters

area_num = 4

# bin number for ra & dec of each exposure
ra_bin_num = 10
dec_bin_num = 10

deg2arc = numpy.pi/180

ra_idx = 29
dec_idx = 30
redshift_idx = 10
redshift_bin_num = 6
red_shift_bin = [0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3]

# star number on each chip
nstar_idx = 21
nstar_thresh = 12

# selection function
flux2_alt_idx = 28
flux2_alt_thresh = 3

# selection on the bad pixels
imax_idx = 22
jmax_idx = 23
imax_thresh = 48
jmax_thresh = 48

# field distortion
gf1_idx = 31
gf2_idx = 32
gf1_thresh = 0.005
gf2_thresh = 0.0075

# shear estimators
mg1_idx = 33
mg2_idx = 34
mn_idx = 35
mu_idx = 36
mv_idx = 37

# exposure label
expo_idx = 38

fourier_cata_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/fourier_cata/hdf5_cata"
result_cata_path = "/mnt/ddnfs/data_users/hkli/CFHT/correlation/cata"

cmd = argv[1]
if cmd == "prepare":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()


    fields, field_name = tool_box.field_dict(fourier_cata_path + "/nname_avail.dat")
    if rank == 0:
        print("Prepare catalog files")
        print(len(field_name))

    field_name_sub = tool_box.alloc(field_name, cpus, method="order")[rank]

    field_avail_sub = []
    field_avail_nm_sub = []


    for fns in field_name_sub:

        field_dst_path = result_cata_path + "/%s.hdf5" % fns
        h5f_dst = h5py.File(field_dst_path,"w")

        # read the the field data
        h5f_src = h5py.File(fourier_cata_path + "/%s/%s.hdf5" % (fns,fns), "r")
        src_data = h5f_src["/field"][()]

        ra = src_data[:,ra_idx]
        dec = src_data[:,dec_idx]

        # find the center and the 4 corners
        ra_min, ra_max = ra.min(), ra.max()
        ra_center = (ra_min + ra_max)/2

        dec_min, dec_max = dec.min(), dec.max()
        dec_center = (dec_min + dec_max)/2

        dec_field = numpy.array([dec_min, dec_min, dec_max, dec_max, dec_center])
        ra_field = numpy.array([ra_min, ra_max, ra_max, ra_min, ra_center])

        block_num = dec_bin_num*ra_bin_num
        dec_bin = numpy.linspace(dec_min, dec_max + 0.01*(dec_max-dec_min), dec_bin_num+1,dtype=numpy.float32)
        ra_bin = numpy.linspace(ra_min, ra_max + 0.01*(ra_max-ra_min), ra_bin_num+1,dtype=numpy.float32)

        block_dec = numpy.zeros((block_num, 4),dtype=numpy.float32)
        block_ra = numpy.zeros((block_num, 4),dtype=numpy.float32)
        for i in range(dec_bin_num):
            for j in range(ra_bin_num):
                ij = i*ra_bin_num + j
                block_dec[ij] = dec_bin[i], dec_bin[i], dec_bin[i+1], dec_bin[i+1]
                block_ra[ij] = ra_bin[j], ra_bin[j+1], ra_bin[j+1], ra_bin[j]

        gal_num_in_block = numpy.zeros((block_num, ), dtype=numpy.intc)
        block_st = numpy.zeros((block_num, ), dtype=numpy.intc)
        block_ed = numpy.zeros((block_num, ), dtype=numpy.intc)

        # read each exposure
        expos_num = h5f_src["/expos_num"][()][0]
        h5f_src.close()

        # selection
        idx1 = src_data[:, nstar_idx] >= nstar_thresh
        idx2 = src_data[:, flux2_alt_idx] >= flux2_alt_thresh
        idx3 = numpy.abs(src_data[:, gf1_idx]) <= gf1_thresh
        idx4 = numpy.abs(src_data[:, gf2_idx]) <= gf2_thresh
        idx5 = src_data[:, imax_idx] <= imax_thresh
        idx6 = src_data[:, jmax_idx] <= jmax_thresh

        idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6

        num_in_z_bin = numpy.zeros((redshift_bin_num,), dtype=numpy.intc)

        for i in range(redshift_bin_num):
            idx_z1 = src_data[:,redshift_idx] >= red_shift_bin[i]
            idx_z2 = src_data[:,redshift_idx] < red_shift_bin[i+1]
            idx_final = idx & idx_z1 & idx_z2

            final_num = idx_final.sum()
            final_data = src_data[idx_final]

            if final_num > 1:
                dst_data = numpy.zeros((final_num, 10), dtype=numpy.float32)

                num_in_z_bin[i] = final_num

                # assign the galaxies to the block in each exposure of each field
                for i in range(dec_bin_num):
                    idx1 = final_data[:, dec_idx] >= dec_bin[i]
                    idx2 = final_data[:, dec_idx] < dec_bin[i+1]
                    idx_d = idx1 & idx2

                    for j in range(ra_bin_num):
                        idx3 = final_data[:, ra_idx] >= ra_bin[j]
                        idx4 = final_data[:, ra_idx] < ra_bin[j+1]

                        idx_block = idx_d & idx3 & idx4

                        ij = i*ra_bin_num+j
                        gal_num_in_block[ij] = idx_block.sum()

                        st = gal_num_in_block[:ij].sum()
                        ed = st + gal_num_in_block[ij]
                        block_st[ij] = st
                        block_ed[ij] = ed

                        dst_data[st:ed, 0] = final_data[:, mg1_idx][idx_block]
                        dst_data[st:ed, 1] = final_data[:, mg2_idx][idx_block]
                        dst_data[st:ed, 2] = final_data[:, mn_idx][idx_block]
                        dst_data[st:ed, 3] = final_data[:, mu_idx][idx_block]
                        dst_data[st:ed, 4] = final_data[:, mv_idx][idx_block]

                        dst_data[st:ed, 5] = final_data[:, dec_idx][idx_block]
                        dst_data[st:ed, 7] = final_data[:, ra_idx][idx_block]
                        dst_data[st:ed, 9] = final_data[:, redshift_idx][idx_block]
                        dst_data[st:ed, 8] = final_data[:, expo_idx][idx_block]

                dst_data[:, 6] = numpy.cos(final_data[:, dec_idx]*deg2arc)
                h5f_dst["/field"] = dst_data
        h5f_dst["/num_in_zbin"] = num_in_z_bin
        



        # target_label = 0
        # source_num_in_expo = [[], []]
        #
        # for iexpo in range(expos_num):
        #
        #     src_expo_data = h5f_src["/expo_%d"%iexpo][()]
        #
        #     src_expo_name = int(h5f_src["/expo_%d"%iexpo].attrs["exposure_name"].split("p")[0])
        #     src_num = src_expo_data.shape[0]
        #
        #     # selection
        #     idx1 = src_expo_data[:, nstar_idx] >= nstar_thresh
        #     idx2 = src_expo_data[:, flux2_alt_idx] >= flux2_alt_thresh
        #     idx3 = numpy.abs(src_expo_data[:, gf1_idx]) <= gf1_thresh
        #     idx4 = numpy.abs(src_expo_data[:, gf2_idx]) <= gf2_thresh
        #     idx5 = src_expo_data[:,imax_idx] <= imax_thresh
        #     idx6 = src_expo_data[:,jmax_idx] <= jmax_thresh
        #
        #     idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6
        #     select_num = idx.sum()
        #
        #     if select_num > 1:
        #
        #         dec_s = src_expo_data[:, dec_idx][idx]
        #         ra_s = src_expo_data[:,ra_idx][idx]
        #         redshift_s = src_expo_data[:,redshift_idx][idx]
        #
        #         # shear estimators
        #         mg1_s = src_expo_data[:,mg1_idx][idx]
        #         mg2_s = src_expo_data[:,mg2_idx][idx]
        #         mn_s = src_expo_data[:,mn_idx][idx]
        #         mu_s = src_expo_data[:,mu_idx][idx]
        #         mv_s = src_expo_data[:,mv_idx][idx]
        #
        #         dst_data = numpy.zeros((select_num, 9), dtype=numpy.float32)
        #         expo_label = numpy.zeros((select_num, 1), dtype=numpy.intc) + src_expo_name
        #
        #         source_num_in_expo[0].append(select_num)
        #         source_num_in_expo[1].append(src_expo_name)
        #
        #         # assign the galaxies to the block in each exposure of each field
        #         for i in range(dec_bin_num):
        #             idx1 = dec_s >= dec_bin[i]
        #             idx2 = dec_s < dec_bin[i+1]
        #             idx_d = idx1 & idx2
        #
        #             for j in range(ra_bin_num):
        #                 idx3 = ra_s >= ra_bin[j]
        #                 idx4 = ra_s < ra_bin[j+1]
        #
        #                 idx_block = idx_d & idx3 & idx4
        #
        #                 ij = i*ra_bin_num+j
        #                 gal_num_in_block[ij] = idx_block.sum()
        #
        #                 st = gal_num_in_block[:ij].sum()
        #                 ed = st + gal_num_in_block[ij]
        #                 block_st[ij] = st
        #                 block_ed[ij] = ed
        #
        #                 dst_data[st:ed, 0] = mg1_s[idx_block]
        #                 dst_data[st:ed, 1] = mg2_s[idx_block]
        #                 dst_data[st:ed, 2] = mn_s[idx_block]
        #                 dst_data[st:ed, 3] = mu_s[idx_block]
        #                 dst_data[st:ed, 4] = mv_s[idx_block]
        #
        #                 dst_data[st:ed, 5] = dec_s[idx_block]
        #                 dst_data[st:ed, 6] = ra_s[idx_block]
        #                 dst_data[st:ed, 7] = redshift_s[idx_block]
        #
        #         dst_data[:, 8] = numpy.cos(dst_data[:,5]*deg2arc)
        #         dst_data[:, 9] = numpy.cos(dst_data[:,5]*deg2arc)
        #
        #         h5f_dst["/expo_%d/mg1"%target_label] = dst_data[:, 0]
        #         h5f_dst["/expo_%d/mg2"%target_label] = dst_data[:, 1]
        #         h5f_dst["/expo_%d/mn"%target_label] = dst_data[:, 2]
        #         h5f_dst["/expo_%d/mu"%target_label] = dst_data[:, 3]
        #         h5f_dst["/expo_%d/mv"%target_label] = dst_data[:, 4]
        #
        #         h5f_dst["/expo_%d/dec"%target_label] = dst_data[:,5]
        #         h5f_dst["/expo_%d/ra"%target_label] = dst_data[:,6]
        #         h5f_dst["/expo_%d/redshift"%target_label] = dst_data[:,7]
        #         h5f_dst["/expo_%d/cos_dec" % target_label] = dst_data[:, 8]
        #
        #         h5f_dst["/expo_%d/gal_num_in_block"%target_label] = gal_num_in_block
        #         h5f_dst["/expo_%d/block_start"%target_label] = block_st
        #         h5f_dst["/expo_%d/block_end"%target_label] = block_ed
        #         h5f_dst["/expo_%d/block_dec"%target_label] = block_dec
        #         h5f_dst["/expo_%d/block_ra"%target_label] = block_ra
        #
        #         h5f_dst["/expo_%d/ra_bin"%target_label] = ra_bin
        #         h5f_dst["/expo_%d/dec_bin"%target_label] = dec_bin
        #
        #         h5f_dst["/expo_%d"%target_label].attrs["exposure_name"] = src_expo_name
        #
        #         if target_label == 0:
        #             stack_data = dst_data
        #             stack_expo_label = expo_label
        #         else:
        #             stack_data = numpy.row_stack((stack_data, dst_data))
        #             stack_expo_label = numpy.row_stack((stack_expo_label, expo_label))
        #
        #         target_label += 1
        #         # print(rank,"%s-%s   %d    %d"%(field_dst_path,src_expo_name, select_num, src_num))
        # # start & end label of each exposure in the stacked data (field)
        # expo_st = []
        # expo_ed = []
        # for iexpo in range(target_label):
        #     expo_st.append(sum(source_num_in_expo[1][:iexpo]))
        #     expo_ed.append(expo_st[iexpo] + source_num_in_expo[1][iexpo])
        #
        # h5f_src.clear()
        #
        # if target_label < 1:
        #     h5f_dst.close()
        #     os.remove(field_dst_path)
        # else:
        #     field_avail_sub.append("%s\t%d\t%d\n"%(field_dst_path, target_label, expos_num))
        #     field_avail_nm_sub.append("%s\n"%fns)
        #
        #     h5f_dst["/expo_start"] = numpy.array(expo_st,dtype=numpy.intc)
        #     h5f_dst["/expo_end"] = numpy.array(expo_ed,dtype=numpy.intc)
        #     h5f_dst["/field_ra"] = ra_field
        #     h5f_dst["/field_dec"] = dec_field
        #
        #     h5f_dst["/field_expo_label"] = stack_expo_label
        #     h5f_dst["/total_gal_num"] = numpy.array([stack_data.shape[0]],dtype=numpy.intc)
        #     h5f_dst["/field"] = stack_data
        #     h5f_dst["/expos_num"] = numpy.array([target_label])
        #
        #     h5f_dst["/expos_info"] = numpy.array(source_num_in_expo, dtype=numpy.intc)
        #     h5f_dst.close()

    comm.Barrier()

    field_avail_list = comm.gather(field_avail_sub, root=0)
    field_avail_nm_list = comm.gather(field_avail_nm_sub, root=0)

    if rank == 0:
        buffer_field = []
        for fns in field_avail_list:
            buffer_field.extend(fns)

        buffer_field_nm = []
        for fns in field_avail_nm_list:
            buffer_field_nm.extend(fns)

        field_name_buffer = []
        field_dst_path_buffer = []

        for i in range(1, area_num+1):
            field_dst_path_buffer_sub = []
            for fns in buffer_field:
                if "w%d"%i in fns:
                    field_dst_path_buffer.append(fns)
                    field_dst_path_buffer_sub.append(fns)
            with open(result_cata_path + "/source_list_%d.dat"%(i-1),"w") as f:
                f.writelines(field_dst_path_buffer_sub)
            print(len(field_dst_path_buffer_sub))

            for fns in buffer_field_nm:
                if "w%d"%i in fns:
                    field_name_buffer.append(fns)

        with open(result_cata_path + "/fields.dat", "w") as f:
            f.writelines(field_name_buffer)

        with open(result_cata_path + "/source_list_all.dat", "w") as f:
            f.writelines(field_dst_path_buffer)


    comm.Barrier()


elif cmd == "stack":
    # stack the fields

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    result_path = result_cata_path + "/total.hdf5"

    if rank == 0:
        print("Stack field files")
        h5f = h5py.File(result_path, "w")
        h5f.close()

    for i in range(area_num):

        with open(result_cata_path + "/source_list_%d.dat" % i) as f:
            fields = f.readlines()
        sub_area_list = []
        for nm in fields:
            sub_area_list.append(nm.split("\t")[0])
        # print(sub_area_list)
        my_sub_area_list = tool_box.alloc(sub_area_list, cpus, method="order")[rank]
        # print(rank, i, len(my_sub_area_list))
        gal_num = []
        if len(my_sub_area_list) > 0:
            for tag, fns in enumerate(my_sub_area_list):

                h5f = h5py.File(fns, "r")

                temp1 = h5f["/field"][()][:, 5:7]
                temp2 = h5f["/field_expo_label"][()]

                num = h5f["/total_gal_num"][()][0]
                expo_info = h5f["/expos_info"][()]
                # print(temp1.shape)
                if num != temp1.shape[0] or expo_info[0].sum() != num:
                    print("Wrong %s" % fns, num,temp1.shape[0], expo_info[0].sum())

                if tag == 0:
                    data1 = temp1
                    data2 = temp2
                else:
                    data1 = numpy.row_stack((data1, temp1))
                    data2 = numpy.row_stack((data2, temp2))
                h5f.close()

                gal_num.append(num)
            sp1 = data1.shape
            sp2 = data2.shape
        else:
            sp1 = (0, 0)
            sp2 = (0, 0)
        # print(sp1, sp2)

        sp1_total = comm.gather(sp1, root=0)
        sp2_total = comm.gather(sp2, root=0)

        # print(i,rank, data.shape, data.dtype, sp, data[0,:5])
        comm.Barrier()

        if rank > 0 and sp1[0] > 0:
            comm.Send([data1, MPI.FLOAT], dest=0, tag=rank)
        else:
            # print(sp1_total)
            for ir in range(1, cpus):
                if sp1_total[ir][0] > 0:
                    recv_buf = numpy.empty(sp1_total[ir], dtype=numpy.float32)
                    comm.Recv(recv_buf, source=ir, tag=ir)
                    data1 = numpy.row_stack((data1, recv_buf))

            if i == 0:
                stack_data1 = data1
            else:
                stack_data1 = numpy.row_stack((stack_data1, data1))

            h5f = h5py.File(result_path, "r+")
            h5f["/w%d" % i] = data1
            h5f["/w%d" % i].attrs["info"] = ["dec ra"]
            h5f.close()
        comm.Barrier()

        if rank > 0 and sp2[0] > 0:
            comm.Send([data2, MPI.INT], dest=0, tag=rank)
        else:
            for ir in range(1, cpus):
                if sp2_total[ir][0] > 0:
                    recv_buf = numpy.empty(sp2_total[ir], dtype=numpy.intc)
                    comm.Recv(recv_buf, source=ir, tag=ir)
                    data2 = numpy.row_stack((data2, recv_buf))

            if i == 0:
                stack_data2 = data2
            else:
                stack_data2 = numpy.row_stack((stack_data2, data2))

            h5f = h5py.File(result_path, "r+")
            h5f["/w%d_expo_label" % i] = data2
            h5f["/w%d_expo_label" % i].attrs["info"] = ["exposure labels"]
            h5f.close()
        comm.Barrier()

        num_buffer = comm.gather(gal_num, root=0)
        if rank == 0:
            source_num_in_field = []
            for snb in num_buffer:
                if len(snb) > 0:
                    source_num_in_field.extend(snb)

            fnum = len(source_num_in_field)
            source_num_in_field = numpy.array(source_num_in_field, dtype=numpy.intc).reshape(fnum, 1)
            if source_num_in_field.sum() != data1.shape[0]:
                print("Wrong-total_num !")

            if i == 0:
                stack_source_num_in_field = source_num_in_field
            else:
                stack_source_num_in_field = numpy.row_stack((stack_source_num_in_field, source_num_in_field))

            h5f = h5py.File(result_path, "r+")
            h5f["/w%d_gal_num" % i] = source_num_in_field
            h5f["/w%d_gal_num" % i].attrs["info"] = ["gal num of each field"]
            h5f.close()

        comm.Barrier()

    if rank == 0:
        h5f = h5py.File(result_path, "r+")

        h5f["/total/data"] = stack_data1
        h5f["/total/data"].attrs["info"] = ["dec ra"]

        h5f["/total/expo_label"] = stack_data2
        h5f["/total/expo_label"].attrs["info"] = ["exposure labels"]

        h5f["/total/gal_num"] = stack_source_num_in_field
        h5f["/total/gal_num"].attrs["info"] = ["gal num of each field"]

        h5f.close()
    comm.Barrier()

elif cmd == "segment":
    # Kmeans method for classification for jackknife

    t1 = time.time()

    njobs = int(argv[2])
    ncent = 200

    h5f = h5py.File(result_cata_path + "/total.hdf5", "r")

    total = h5f["/total/data"][()]
    expo_label = h5f["/total/expo_label"][()]
    gal_num = h5f["/total/gal_num"][()]
    h5f.close()

    y_pred = KMeans(n_clusters=ncent, random_state=numpy.random.randint(1, 100000), n_jobs=njobs).fit_predict(total)

    h5f = h5py.File(result_cata_path + "/total.hdf5", "r+")
    h5f["/total/area_labels_mini"] = y_pred
    h5f["/total/area_labels_mini"].attrs["info"] = ["labels from Kmeans for jackknife"]
    h5f.close()

    t2 = time.time()
    print(t2-t1)