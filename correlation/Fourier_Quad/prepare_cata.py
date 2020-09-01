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
# theta bin
theta_bin_num = 5
theta_bin = tool_box.set_bin_log(0.8, 60, theta_bin_num+1).astype(numpy.float32)

# bin number for ra & dec of each exposure
deg2arcmin = 60
deg2rad = numpy.pi/180

grid_size = 15 #arcmin

ra_idx = 0
dec_idx = 1
redshift_idx = 10
redshift_bin_num = 6
redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)

# chi guess bin for PDF_SYM
chi_guess_num = 32
chi_guess_bin = tool_box.set_bin_log(10**(-7), 10**(-3), chi_guess_num).astype(numpy.float32)
cor_gg_len = 524288
mg_bin_num = 10

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

#
# fourier_cata_path = "/coma/hklee/fourier_cata/original_cata"
# result_cata_path = "/coma/hklee/correlation/cata"
fourier_cata_path = "/mnt/perc/hklee/CFHT/catalog/fourier_cata/original_cata"
result_cata_path = "/mnt/perc/hklee/CFHT/correlation/cata"

cmd = argv[1]

if cmd == "correlation":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()
    # rank = 0
    # cpus = 1

    # for correlation calculation
    if rank == 0:

        field_name = []
        with open(fourier_cata_path + "/nname_avail.dat", "r") as f:
            f_lines = f.readlines()
        for ff in f_lines:
            field_name.append(ff.split("\n")[0])

        # set up bins for G1(or G2) for PDF_SYM
        for i in range(10):
            h5f_src = h5py.File(fourier_cata_path + "/%s/result/%s.hdf5" % (field_name[i], field_name[i]), "r")
            temp = h5f_src["/field"][()][:, mg1_idx:mg2_idx + 1]
            if i == 0:
                src_data = temp
            else:
                src_data = numpy.row_stack((src_data, temp))

        mg_bin = tool_box.set_bin(src_data[:, 0], mg_bin_num, 100000)


        h5f_cor = h5py.File(result_cata_path + "/gg_cor.hdf5", "w")

        for i in range(chi_guess_num):

            mean = [0, 0]

            cov = [[chi_guess_bin[i] * 2, chi_guess_bin[i]],
                   [chi_guess_bin[i], chi_guess_bin[i] * 2]]

            gg = tool_box.rand_gauss2n(cor_gg_len, mean, cov).astype(dtype=numpy.float32)

            if i == 0:
                gg_1 = gg[0]
                gg_2 = gg[1]
            else:
                gg_1 = numpy.column_stack((gg_1, gg[0]))
                gg_2 = numpy.column_stack((gg_2, gg[1]))
        print(gg_1.shape)
        h5f_cor["/g11"] = gg_1.reshape((chi_guess_num * cor_gg_len,))
        h5f_cor["/g22"] = gg_2.reshape((chi_guess_num * cor_gg_len,))
        h5f_cor["/chi_guess"] = chi_guess_bin
        h5f_cor["/theta_bin"] = theta_bin
        h5f_cor["/theta_bin_num"] = numpy.array([theta_bin_num], dtype=numpy.intc)
        h5f_cor["/redshift_bin"] = redshift_bin
        h5f_cor["/redshift_bin_num"] = numpy.array([redshift_bin_num], dtype=numpy.intc)
        h5f_cor["/mg_bin"] = mg_bin.astype(dtype=numpy.float32)

        h5f_cor.close()

    comm.Barrier()


if cmd == "prepare_exposure_wise":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    field_name = []
    with open(fourier_cata_path + "/nname_avail.dat","r") as f:
        f_lines = f.readlines()
    for ff in f_lines:
        field_name.append(ff.split("\n")[0])

    fields, field_name = tool_box.field_dict(fourier_cata_path + "/nname_all.dat")

    total_expos = []
    for fnm in field_name:
        if "w1" in fnm:
            sub_exp = fields[fnm]
            for expo in sub_exp:
                total_expos.append((fnm, expo))



    my_expos = tool_box.alloc(total_expos,cpus, method="order")[rank]
    expo_avail_sub = []

    for expo in my_expos:
        h5f_src = h5py.File(fourier_cata_path + "/%s/result/%s_all.hdf5"%(expo[0], expo[1]),"r")
        data = h5f_src["/data"][()]
        h5f_src.close()

        # selection
        idx1 = data[:, nstar_idx] >= nstar_thresh
        idx2 = data[:, flux2_alt_idx] >= flux2_alt_thresh
        idx3 = numpy.abs(data[:, gf1_idx]) <= gf1_thresh
        idx4 = numpy.abs(data[:, gf2_idx]) < gf2_thresh
        idx5 = data[:, imax_idx] <= imax_thresh
        idx6 = data[:, jmax_idx] <= jmax_thresh
        idx_7 = data[:, redshift_idx] >= redshift_bin[0]
        idx_8 = data[:, redshift_idx] < redshift_bin[redshift_bin_num]

        idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx_7 & idx_8

        src_num = idx.sum()

        if src_num > 0:
            src_data = data[idx]

            src_data[:, ra_idx] = src_data[:, ra_idx] * deg2arcmin
            src_data[:, dec_idx] = src_data[:, dec_idx] * deg2arcmin
            ra = src_data[:, ra_idx]
            dec = src_data[:, dec_idx]

            # find the center and the 4 corners
            ra_min, ra_max = ra.min(), ra.max()
            ra_center = (ra_min + ra_max) / 2
            dra = (ra_max - ra_min) / 2
            dec_min, dec_max = dec.min(), dec.max()
            dec_center = (dec_min + dec_max) / 2
            cos_dec_center = numpy.cos(dec_center / deg2arcmin * deg2rad)
            ddec = (dec_max - dec_min) / 2

            expo_pos = numpy.array([ra_center, dec_center, dra, ddec,
                                     numpy.sqrt((dra * cos_dec_center) ** 2 + ddec ** 2), cos_dec_center],
                                    dtype=numpy.float32)

            dst_data = numpy.zeros((src_num, 9), dtype=numpy.float32)

            redshift_label = numpy.zeros((src_num,), dtype=numpy.intc)

            test_mask = numpy.zeros((src_num, ), dtype=numpy.intc)

            total_num_in_z_bin = numpy.zeros((redshift_bin_num, ), dtype=numpy.intc)
            iz_st = numpy.zeros((redshift_bin_num, ), dtype=numpy.intc)
            iz_ed = numpy.zeros((redshift_bin_num, ), dtype=numpy.intc)

            for iz in range(redshift_bin_num):
                idx_z1 = src_data[:,redshift_idx] >= redshift_bin[iz]
                idx_z2 = src_data[:,redshift_idx] < redshift_bin[iz+1]
                idx = idx_z1 & idx_z2

                total_num_in_z_bin[iz] = idx.sum()

                st = total_num_in_z_bin[:iz].sum()
                ed = st + total_num_in_z_bin[iz]
                # print(iz, st, ed)

                iz_st[iz] = st
                iz_ed[iz] = ed

                dst_data[st:ed, 0] = src_data[:, mg1_idx][idx]
                dst_data[st:ed, 1] = src_data[:, mg2_idx][idx]
                dst_data[st:ed, 2] = src_data[:, mn_idx][idx]
                dst_data[st:ed, 3] = src_data[:, mu_idx][idx]
                dst_data[st:ed, 4] = src_data[:, mv_idx][idx]

                dst_data[st:ed, 5] = src_data[:, ra_idx][idx]
                dst_data[st:ed, 6] = src_data[:, dec_idx][idx]
                dst_data[st:ed, 8] = src_data[:, redshift_idx][idx]
                redshift_label[st:ed] = iz
                test_mask[st:ed] = 1

            if test_mask.sum() != src_num:
                print("%d  %s-%s wrong in mask(%d!=%d)" % (rank, expo[0], expo[1], test_mask.sum(), src_num))
            dst_data[:, 7] = numpy.cos(dst_data[:, 6] / deg2arcmin * deg2rad)
            # print(total_num_in_z_bin)
            expo_dst_path = result_cata_path + "/%s-%s.hdf5" %(expo[0], expo[1])
            h5f_dst = h5py.File(expo_dst_path,"w")

            h5f_dst["/pos"] = expo_pos
            h5f_dst["/redshift_label"] = redshift_label
            h5f_dst["/data"] = dst_data
            h5f_dst["/num_in_zbin"] = total_num_in_z_bin
            h5f_dst["/zbin_st"] = iz_st
            h5f_dst["/zbin_ed"] = iz_ed
            h5f_dst.close()

            expo_avail_sub.append("%s\t%s-%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n"
                                   %(expo_dst_path, expo[0], expo[1], src_num, expo_pos[0],expo_pos[1], expo_pos[2],
                                     expo_pos[3], expo_pos[4], expo_pos[5]))
        # exit()
    comm.Barrier()

    expo_avail_list = comm.gather(expo_avail_sub, root=0)

    if rank == 0:
        buffer_expo = []
        for fns in expo_avail_list:
            buffer_expo.extend(fns)

        with open(result_cata_path + "/source_list.dat", "w") as f:
            f.writelines(buffer_expo)
        # print(buffer_expo)
        print("Prepare catalog files")
        # print(total_expos)
        print(len(field_name), "fields", len(total_expos), " exposures")
    comm.Barrier()


if cmd == "prepare_field_wise":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    fields, field_name = tool_box.field_dict(fourier_cata_path + "/nname_avail.dat")
    if rank == 0:
        print("Prepare catalog files")
        print(len(field_name))

    field_name_sub = tool_box.alloc(field_name, cpus, method="order")[rank]

    field_avail_sub = []

    for fns in field_name_sub:

        field_dst_path = result_cata_path + "/%s.hdf5" % fns
        h5f_dst = h5py.File(field_dst_path,"w")

        # read the the field data
        h5f_src = h5py.File(fourier_cata_path + "/%s/result/%s.hdf5" % (fns,fns), "r")
        src_data = h5f_src["/field"][()]

        src_data[:, ra_idx] = src_data[:,ra_idx]*deg2arcmin
        src_data[:, dec_idx] = src_data[:,dec_idx]*deg2arcmin
        ra = src_data[:,ra_idx]
        dec = src_data[:,dec_idx]

        # find the center and the 4 corners
        ra_min, ra_max = ra.min(), ra.max()
        ra_center = (ra_min + ra_max)/2
        dra = (ra_max - ra_min)/2
        dec_min, dec_max = dec.min(), dec.max()
        dec_center = (dec_min + dec_max)/2
        cos_dec_center = numpy.cos(dec_center/deg2arcmin*deg2rad)
        ddec = (dec_max - dec_min)/2

        field_pos = numpy.array([ra_center, dec_center, dra, ddec,
                                 numpy.sqrt((dra*cos_dec_center)**2 + ddec**2),cos_dec_center],dtype=numpy.float32)

        grid_pad = grid_size*0.1
        dec_bin_num = int((dec_max + grid_pad - (dec_min - grid_pad))/grid_size)+1
        ra_bin_num = int((ra_max + grid_pad - (ra_min - grid_pad))/grid_size)+1

        block_num = dec_bin_num*ra_bin_num

        dec_bin = numpy.array([dec_min - grid_pad + grid_size*i for i in range(dec_bin_num+1)], dtype=numpy.float32)
        ra_bin = numpy.array([ra_min - grid_pad + grid_size*i for i in range(ra_bin_num+1)], dtype=numpy.float32)

        # print(dec_bin)
        # print(ra_bin)
        # print(dec_min, dec_max, dec_bin_num)
        # print(ra_min, ra_max, ra_bin_num)
        block_pos = numpy.zeros((block_num, 4),dtype=numpy.float32)
        for i in range(dec_bin_num):
            for j in range(ra_bin_num):
                ij = i*ra_bin_num + j
                block_pos[ij,:2] = (ra_bin[j]+ra_bin[j+1])/2, (dec_bin[i] + dec_bin[i+1])/2
        block_pos[:,2] = numpy.cos(block_pos[:,1]/deg2arcmin*deg2rad)
        block_pos[:,3] = numpy.sqrt(grid_size/2*grid_size/2 + grid_size/2*block_pos[:,2]*grid_size/2*block_pos[:,2])

        h5f_dst["/block_cen_ra"] = block_pos[:,0]
        h5f_dst["/block_cen_dec"] = block_pos[:,1]
        h5f_dst["/block_cen_cos_dec"] = block_pos[:,2]
        h5f_dst["/block_delta_len"] = block_pos[:,3]


        gal_num_in_block = numpy.zeros((block_num, ), dtype=numpy.intc)
        block_st = numpy.zeros((block_num, ), dtype=numpy.intc)
        block_ed = numpy.zeros((block_num, ), dtype=numpy.intc)

        # read each exposure
        expos_num = h5f_src["/expos_num"][()][0]
        h5f_src.close()

        field_avail_sub.append("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n"
                               %(field_dst_path, fns, expos_num, field_pos[0],
                                  field_pos[1], field_pos[2], field_pos[3], field_pos[4], field_pos[5]))

        # selection
        idx1 = src_data[:, nstar_idx] >= nstar_thresh
        idx2 = src_data[:, flux2_alt_idx] >= flux2_alt_thresh
        idx3 = numpy.abs(src_data[:, gf1_idx]) <= gf1_thresh
        idx4 = numpy.abs(src_data[:, gf2_idx]) <= gf2_thresh
        idx5 = src_data[:, imax_idx] <= imax_thresh
        idx6 = src_data[:, jmax_idx] <= jmax_thresh

        idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6

        total_num_in_z_bin = numpy.zeros((redshift_bin_num,), dtype=numpy.intc)

        for iz in range(redshift_bin_num):
            idx_z1 = src_data[:,redshift_idx] >= redshift_bin[iz]
            idx_z2 = src_data[:,redshift_idx] < redshift_bin[iz+1]
            idx_final = idx & idx_z1 & idx_z2

            final_num = idx_final.sum()
            final_data = src_data[idx_final]
            # print(fns, iz, final_num, final_data.shape)
            if final_num > 1:
                dst_data = numpy.zeros((final_num, 11), dtype=numpy.float32)
                test_mask = numpy.zeros((final_num, ), dtype=numpy.intc)

                total_num_in_z_bin[iz] = final_num

                # assign the galaxies to the block in each exposure of each field
                for i in range(dec_bin_num):
                    idx1 = final_data[:, dec_idx] >= dec_bin[i]
                    idx2 = final_data[:, dec_idx] < dec_bin[i+1]
                    idx_d = idx1 & idx2
                    # print(idx_d.sum())
                    for j in range(ra_bin_num):
                        idx3 = final_data[:, ra_idx] >= ra_bin[j]
                        idx4 = final_data[:, ra_idx] < ra_bin[j+1]

                        idx_block = idx_d & idx3 & idx4

                        ij = i*ra_bin_num+j
                        gal_num_in_block[ij] = idx_block.sum()
                        # print(iz, ij, idx_block.sum())
                        st = gal_num_in_block[:ij].sum()
                        ed = st + gal_num_in_block[ij]
                        block_st[ij] = st
                        block_ed[ij] = ed

                        dst_data[st:ed, 0] = final_data[:, mg1_idx][idx_block]
                        dst_data[st:ed, 1] = final_data[:, mg2_idx][idx_block]
                        dst_data[st:ed, 2] = final_data[:, mn_idx][idx_block]
                        dst_data[st:ed, 3] = final_data[:, mu_idx][idx_block]
                        dst_data[st:ed, 4] = final_data[:, mv_idx][idx_block]

                        dst_data[st:ed, 5] = final_data[:, ra_idx][idx_block]
                        dst_data[st:ed, 6] = final_data[:, dec_idx][idx_block]
                        dst_data[st:ed, 8] = final_data[:, redshift_idx][idx_block]

                        dst_data[st:ed, 9] = final_data[:, expo_idx][idx_block] - 1
                        dst_data[st:ed, 10] = ij

                        test_mask[st:ed] = 1
                # print(gal_num_in_block)
                if test_mask.sum() != final_num:
                    print("%d  %s z%d wrong in mask(%d!=%d)"%(rank, fns, iz, test_mask.sum(), final_num))
                dst_data[:, 7] = numpy.cos(dst_data[:, 6]/deg2arcmin * deg2rad)

                h5f_dst["/z%d/field"%iz] = dst_data[:,:9]

                h5f_dst["/z%d/gal_num_in_block"%iz] = gal_num_in_block
                h5f_dst["/z%d/block_st"%iz] = block_st
                h5f_dst["/z%d/block_ed"%iz] = block_ed
                h5f_dst["/z%d/ra_bin"%iz] = ra_bin
                h5f_dst["/z%d/dec_bin"%iz] = dec_bin
                h5f_dst["/z%d/exposure_label"%iz] = (dst_data[:,9] - dst_data[:,9].min()).astype(dtype=numpy.intc)
                h5f_dst["/z%d/block_label"%iz] = dst_data[:,10].astype(dtype=numpy.intc)

        h5f_dst["/field_pos"] = field_pos
        h5f_dst["/expo_num"] = numpy.array([expos_num], dtype=numpy.intc)
        h5f_dst["/total_num_in_zbin"] = total_num_in_z_bin
        h5f_dst.close()


    comm.Barrier()
    #
    field_avail_list = comm.gather(field_avail_sub, root=0)


    if rank == 0:
        buffer_field = []
        for fns in field_avail_list:
            buffer_field.extend(fns)

        with open(result_cata_path + "/source_list.dat", "w") as f:
            f.writelines(buffer_field)

        # set up bins for G1(or G2) for PDF_SYM
        for i in range(5):
            h5f_src = h5py.File(fourier_cata_path + "/%s/result/%s.hdf5" % (field_name[i], field_name[i]), "r")
            temp = h5f_src["/field"][()][:,mg1_idx:mg2_idx+1]
            if i == 0:
                src_data = temp
            else:
                src_data = numpy.row_stack((src_data, temp))

        mg_bin = tool_box.set_bin(temp[:,0], mg_bin_num, 100000)

        h5f_cor = h5py.File(result_cata_path + "/gg_cor.hdf5", "r+")

        h5f_cor["/mg_bin"] = mg_bin.astype(dtype=numpy.float32)
        h5f_cor.close()

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
