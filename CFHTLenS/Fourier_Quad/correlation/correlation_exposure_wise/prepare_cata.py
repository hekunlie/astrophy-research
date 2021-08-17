import os
my_home = os.popen("echo $HK_MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from mpi4py import MPI
import hk_tool_box
import warnings
from sklearn.cluster import KMeans
import time
import hk_FQlib
# import hk_c4py

warnings.filterwarnings('error')

# start
time_start = hk_tool_box.get_time_now()

# parameters

area_num = 4
# theta bin
theta_bin_num = 5
# theta_bin = hk_tool_box.set_bin_log(1, 128, theta_bin_num+1).astype(numpy.float32)
theta_bin = hk_tool_box.set_bin_log(5, 60, theta_bin_num+1).astype(numpy.float32)

# bin number for ra & dec of each exposure
deg2arcmin = 60
deg2rad = numpy.pi/180

grid_size = 15 #arcmin

# position in CFHTLens catalog
ra_idx = 0
dec_idx = 1

redshift_idx = 10

redshift_sep_thresh = 0.01
redshift_bin_num = 6
redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)
# redshift_bin = numpy.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4],dtype=numpy.float32)

# chi guess bin for PDF_SYM
chi_guess_num = 40
inv = [chi_guess_num-1-i for i in range(chi_guess_num)]
chi_guess_bin_p = hk_tool_box.set_bin_log(10**(-8), 10**(-3), chi_guess_num).astype(numpy.float32)
chi_guess_bin = numpy.zeros((2*chi_guess_num, ), dtype=numpy.float32)
chi_guess_bin[:chi_guess_num] = -chi_guess_bin_p[inv]
chi_guess_bin[chi_guess_num:] = chi_guess_bin_p

chi_guess_num = int(chi_guess_num*2)
cor_gg_len = 4#int(chi_guess_num)
cor_gg_len_mid = int(cor_gg_len/2)
mg_bin_num = 10


# star number on each chip
nstar_idx = 21
nstar_thresh = 0

# selection function
flux2_alt_idx = 28
flux2_alt_thresh = 2.5

# selection on the bad pixels
imax_idx = 22
jmax_idx = 23
imax_thresh = 48
jmax_thresh = 48

# the chip labels
ichip_idx = 21

# flux weighted centroid
xc_idx = 23
yc_idx = 24

# field distortion
gf1_idx = 31
gf2_idx = 32
gf1_thresh = 0.005
gf2_thresh = 0.005

# shear estimators
mg1_idx = 33
mg2_idx = 34
mn_idx = 35
mu_idx = 36
mv_idx = 37

# about PhotoZ
odd_idx = 40
odd_thresh = 0.0
redshift_e_idx = 41
# pz_sum_idx = 42
# pz_sum_thresh = 0.001

# 48/2*0.187/60 arcmin
sep_ang = 0.0748
sep_pix = 12 # pixels
sep_z = 0.2


fourier_cata_path = "/home/hklee/work/CFHT/CFHT_shear_cat_6_14_2021_smoothed"
source_list_path = fourier_cata_path + "/cat_inform/field_avail.dat"

result_cata_path = "/home/hklee/work/CFHT/correlation/smooth_cata/cata_odds_0.0/stack_expo"

cmd_avail = "cmd: correlation, prepare, stack, kmeans, segment\n"

cmd = argv[1]

if cmd == "correlation":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()
    if rank == 0:
        print(cmd_avail)

        print(time_start, "correlation")
        print(source_list_path)
    # rank = 0
    # cpus = 1

    # for correlation calculation
    if rank == 0:

        field_name = []
        with open(source_list_path, "r") as f:
            f_lines = f.readlines()
        for ff in f_lines:
            field_name.append(ff.split("\n")[0])

        # set up bins for G1(or G2) for PDF_SYM
        for i in range(20):
            h5f_src = h5py.File(field_name[i], "r")
            temp = h5f_src["/data"][()][:, mg1_idx:mg2_idx + 1]
            if i == 0:
                src_data = temp
            else:
                src_data = numpy.row_stack((src_data, temp))

        mg_bin = hk_tool_box.set_bin(src_data[:, 0], mg_bin_num, 1000000)

        if not os.path.exists(result_cata_path + "/kmeans"):
            os.mkdir(result_cata_path + "/kmeans")
        h5f_cor = h5py.File(result_cata_path + "/kmeans/gg_cor.hdf5", "w")


        gg_1 = numpy.zeros((cor_gg_len,chi_guess_num),dtype=numpy.float32)
        gg_2 = numpy.zeros_like(gg_1,dtype=numpy.float32)

        #
        # for i in range(chi_guess_num):
        #
        #     mean = [0, 0]
        #
        #     cov = [[numpy.abs(chi_guess_bin[i] * 2), chi_guess_bin[i]],
        #            [chi_guess_bin[i], numpy.abs(chi_guess_bin[i] * 2)]]
        #
        #     gg = hk_tool_box.rand_gauss2n(cor_gg_len, mean, cov).astype(dtype=numpy.float32)
        #
        #     gg_1[:,i] = gg[0]
        #     gg_2[:,i] = gg[1]

        # proposed by Zhenjie Liu
        for i in range(chi_guess_num):
            chi_i = chi_guess_bin[i]
            gg_val = numpy.sqrt(numpy.abs(chi_i))

            gg_1[:cor_gg_len_mid, i] = gg_val
            gg_1[cor_gg_len_mid:, i] = -gg_val

            if chi_i > 0:
                gg_2[:cor_gg_len_mid, i] = gg_val
                gg_2[cor_gg_len_mid:, i] = -gg_val
            else:
                gg_2[:cor_gg_len_mid, i] = -gg_val
                gg_2[cor_gg_len_mid:, i] = gg_val

        h5f_cor["/g11_ori"] = gg_1
        h5f_cor["/g22_ori"] = gg_2
        print(gg_1.shape)

        h5f_cor["/g11"] = gg_1.flatten()
        h5f_cor["/g22"] = gg_2.flatten()

        h5f_cor["/chi_guess"] = chi_guess_bin
        h5f_cor["/theta_bin"] = theta_bin
        h5f_cor["/theta_bin_num"] = numpy.array([theta_bin_num], dtype=numpy.intc)
        h5f_cor["/redshift_bin"] = redshift_bin
        h5f_cor["/redshift_bin_num"] = numpy.array([redshift_bin_num], dtype=numpy.intc)
        h5f_cor["/mg_bin"] = mg_bin.astype(dtype=numpy.float32)

        h5f_cor.close()

    comm.Barrier()


if cmd == "prepare":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    if rank == 0:
        print(cmd_avail)

        print(time_start, "prepare")
        print(source_list_path)

    total_expos = []
    with open(source_list_path, "r") as f:
        f_lines = f.readlines()
    for ff in f_lines:
        total_expos.append(ff.split("\n")[0])

    my_expos = hk_tool_box.alloc(total_expos,cpus, method="order")[rank]
    expo_avail_sub = []
    exception_sub = []

    for cat_path in my_expos:
        expo_name = cat_path.split("/")[-1]

        h5f_src = h5py.File(cat_path,"r")
        data = h5f_src["/data"][()]
        h5f_src.close()

        # labels = c4py.deblend_i(data[:,xc_idx], data[:,yc_idx], data[:,redshift_idx], data[:,ichip_idx], sep_pix, sep_z)

        # selection
        idx1 = data[:, nstar_idx] >= nstar_thresh
        idx2 = data[:, flux2_alt_idx] >= flux2_alt_thresh
        idx3 = numpy.abs(data[:, gf1_idx]) <= gf1_thresh
        idx4 = numpy.abs(data[:, gf2_idx]) <= gf2_thresh
        idx5 = data[:, imax_idx] < imax_thresh
        idx6 = data[:, jmax_idx] < jmax_thresh
        idx_7 = data[:, redshift_idx] >= redshift_bin[0]
        idx_8 = data[:, redshift_idx] < redshift_bin[redshift_bin_num]
        idx_9 = data[:,odd_idx] >= odd_thresh
        # idx_10 = labels > 0
        idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx_7 & idx_8 & idx_9# & idx_10

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
            # G1, G2, N, U, V, RA, DEC, COS(DEC), REDSHIFT, REDSHIFT_E
            dst_data = numpy.zeros((src_num, 10), dtype=numpy.float32)

            dst_data[:, 0] = src_data[:, mg1_idx]
            dst_data[:, 1] = src_data[:, mg2_idx]
            dst_data[:, 2] = src_data[:, mn_idx]
            # the signs of U & V are different with that in the paper
            dst_data[:, 3] = -src_data[:, mu_idx]
            dst_data[:, 4] = -src_data[:, mv_idx]

            dst_data[:, 5] = src_data[:, ra_idx]
            dst_data[:, 6] = src_data[:, dec_idx]
            dst_data[:, 7] = numpy.cos(dst_data[:, 6] / deg2arcmin * deg2rad)
            dst_data[:, 8] = src_data[:, redshift_idx]
            dst_data[:, 9] = src_data[:, redshift_e_idx]

            redshift_label = numpy.zeros((src_num,), dtype=numpy.intc)

            test_mask = numpy.zeros((src_num, ), dtype=numpy.intc)

            for iz in range(redshift_bin_num):
                idx_z1 = src_data[:,redshift_idx] >= redshift_bin[iz]
                idx_z2 = src_data[:,redshift_idx] < redshift_bin[iz+1]
                idx = idx_z1 & idx_z2

                redshift_label[idx] = iz
                test_mask[idx] = 1

            if test_mask.sum() != src_num:
                log_inform = "%d %s wrong in mask(%d!=%d)\n" % (rank, expo_name, test_mask.sum(), src_num)
                exception_sub.append(log_inform)

            expo_dst_path = result_cata_path + "/%s" %expo_name
            h5f_dst = h5py.File(expo_dst_path,"w")
            h5f_dst["/pos"] = expo_pos
            h5f_dst["/redshift_label"] = redshift_label
            h5f_dst["/data"] = dst_data
            h5f_dst.close()

            expo_avail_sub.append("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n"
                                   %(expo_dst_path, expo_name.split(".")[0], src_num, expo_pos[0],expo_pos[1], expo_pos[2],
                                     expo_pos[3], expo_pos[4], expo_pos[5]))
        # exit()
    comm.Barrier()

    expo_avail_list = comm.gather(expo_avail_sub, root=0)
    exception_collection = comm.gather(exception_sub, root=0)

    if rank == 0:
        buffer_expo = []
        for fns in expo_avail_list:
            buffer_expo.extend(fns)

        with open(result_cata_path + "/source_list.dat", "w") as f:
            f.writelines(buffer_expo)

        exception_all = []
        for ec in exception_collection:
            exception_all.extend(ec)
        log_inform = "%d (%d) exposures\n"%(len(buffer_expo), len(total_expos))
        exception_all.append(log_inform)
        with open("log.dat", "w") as f:
            f.writelines(exception_all)
        print(log_inform)

        test_nm1 = []
        test_nm2 = []
        for ff in total_expos:
            test_nm1.append(ff.split("/")[-1])
        for ff in buffer_expo:
            test_nm2.append(ff.split()[1]+".hdf5")
        for t1 in test_nm1:
            if t1 not in test_nm2:
                print(t1)
                # pass
        # print(test_nm1)
        # print(test_nm2)
    comm.Barrier()


elif cmd == "stack":
    # collection

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    if rank == 0:
        print(cmd_avail)

        print(time_start, "stack")

    expos = []
    expo_name = []
    gal_num = []
    with open(result_cata_path + "/source_list.dat", "r") as f:
        conts = f.readlines()
    for nm in conts:
        informs = nm.split()
        expos.append(informs[0])
        expo_name.append(informs[1])
        gal_num.append(int(informs[2]))
    expos_num = len(expos)
    if rank == 0:
        print(expos_num, " exposures")

    expos_labels = [i for i in range(expos_num)]
    my_sub_area_list = hk_tool_box.alloc(expos, cpus)[rank]
    my_sub_expos_labels = hk_tool_box.alloc(expos_labels, cpus)[rank]
    my_sub_gal_num = hk_tool_box.alloc(gal_num, cpus)[rank]
    # print(rank, i, len(my_sub_area_list))

    if len(my_sub_area_list) > 0:
        for tag, expo_path in enumerate(my_sub_area_list):

            h5f = h5py.File(expo_path, "r")
            temp = h5f["/data"][()]
            h5f.close()

            temp_label = numpy.zeros((my_sub_gal_num[tag],1),dtype=numpy.intc)
            temp_label[:,0] = my_sub_expos_labels[tag]

            # Nan check
            idx = numpy.isnan(temp)
            if idx.sum() > 0:
                print("Find Nan ", expo_path)
            if tag == 0:
                stack_data = temp
                stack_labels = temp_label
            else:
                stack_data = numpy.row_stack((stack_data, temp))
                stack_labels = numpy.row_stack((stack_labels, temp_label))

        sp = stack_data.shape
    else:
        sp = (0, 0)

    sp_total = comm.gather(sp, root=0)
    # print(i,rank, data.shape, data.dtype, sp, data[0,:5])
    comm.Barrier()

    if rank > 0 and sp[0] > 0:
        comm.Send([stack_data, MPI.FLOAT], dest=0, tag=rank)
    else:
        for ir in range(1, cpus):
            if sp_total[ir][0] > 0:
                recv_buf = numpy.empty(sp_total[ir], dtype=numpy.float32)
                comm.Recv(recv_buf, source=ir, tag=ir)
                stack_data = numpy.row_stack((stack_data, recv_buf))

        # h5f = h5py.File(result_cata_path + "/stack_data.hdf5", "w")
        # h5f["/data"] = stack_data
        # h5f.close()
    comm.Barrier()

    if rank > 0 and sp[0] > 0:
        comm.Send([stack_labels, MPI.INT], dest=0, tag=rank)
    else:
        for ir in range(1, cpus):
            if sp_total[ir][0] > 0:
                recv_buf = numpy.empty((sp_total[ir][0],1), dtype=numpy.intc)
                comm.Recv(recv_buf, source=ir, tag=ir)
                stack_labels = numpy.row_stack((stack_labels, recv_buf))

        # separate it into 4 sub-files
        ra_bin = [1000, 5000, 10000, 15000, 30000]
        # ra_bin = [5000, 10000, 15000, 30000]
        for ir in range(area_num):
            idx1 = stack_data[:,5] >= ra_bin[ir]
            idx2 = stack_data[:,5] < ra_bin[ir+1]
            idx = idx1 & idx2

            h5f = h5py.File(result_cata_path + "/stack_data_%d.hdf5"%ir, "w")
            h5f["/data"] = stack_data[idx]
            h5f["/data_len"] = numpy.array([idx.sum()])
            h5f["/expos_label"] = stack_labels[:,0][idx]
            h5f.close()

        # h5f = h5py.File(result_cata_path + "/stack_data.hdf5", "r+")
        # h5f["/expos_label"] = stack_labels[:,0]
        # h5f["/deblend_label"] = stack_labels[:,1]
        # h5f.close()

        print("Totally %d (%d) galaxies" % (stack_labels.shape[0], sum(gal_num)))

    comm.Barrier()


elif cmd == "kmeans":
    # Kmeans method for classification for jackknife
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()


    if rank == 0:
        print(cmd_avail)

        print(time_start, "kmeans")

    t1 = time.time()

    # 200 sub-samples
    # ncent = [75, 34, 56, 35][rank]
    total_cent = 200
    src_num_each_area = numpy.zeros((area_num,), dtype=numpy.intc)

    if rank < area_num:

        h5f = h5py.File(result_cata_path + "/stack_data_%d.hdf5"%rank, "r")
        data = h5f["/data"][()]
        expos_labels = h5f["/expos_label"][()]
        h5f.close()

        for i in range(area_num):
            h5f = h5py.File(result_cata_path + "/stack_data_%d.hdf5" %i, "r")
            sub_data_len = h5f["/data_len"][()]
            h5f.close()
            src_num_each_area[i] = sub_data_len[0]

        ra_dec = data[:,5:7]

        ncents, ratio = hk_tool_box.even_area(src_num_each_area, area_num, total_cent)

        if rank == 0:
            print(ncents)
            print(ratio)

        ncent = ncents[rank]

        group_pred = KMeans(n_clusters=ncent, random_state=numpy.random.randint(1, 100000)).fit_predict(ra_dec)

        h5f = h5py.File(result_cata_path + "/group_predict_%d.hdf5"%rank, "w")
        h5f["/data"] = group_pred
        h5f["/ra_dec"] = ra_dec
        h5f.close()

        t2 = time.time()
        print("%d Total sub-sample: %d (%d ~ %d). Time: %.2f sec."%(rank, ncent, group_pred.min(), group_pred.max(), t2-t1))
    comm.Barrier()


elif cmd == "segment":
    # assign the source into the artificial exposures

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    if rank == 0:
        print(cmd_avail)

        print(time_start, "segment")

    group_cata_path = result_cata_path + "/kmeans"
    if rank == 0:
        if not os.path.exists(group_cata_path):
            os.makedirs(group_cata_path)
        h5f = h5py.File(result_cata_path + "/bias_test.hdf5", "w")
        h5f.close()
    comm.Barrier()

    # the width of each exposure, arcmin
    expos_field_width = 40
    dra = expos_field_width/2
    ddec = expos_field_width/2

    ncent = []

    expos_avail_sub = []
    expos_count = 0
    for sub_area in range(area_num):
        # read the group labels from Kmeans
        h5f = h5py.File(result_cata_path + "/group_predict_%d.hdf5"%sub_area,"r")
        group_label = h5f["/data"][()]
        h5f.close()

        # the catalog
        h5f = h5py.File(result_cata_path + "/stack_data_%d.hdf5"%sub_area, "r")
        area_data = h5f["/data"][()]
        expos_labels = h5f["/expos_label"][()]
        h5f.close()

        group_num = group_label.max() + 1
        ncent_before = sum(ncent)
        ncent.append(group_num)

        group_list = [i for i in range(group_num)]

        sub_group_list = hk_tool_box.alloc(group_list, cpus)[rank]

        add_bias = numpy.zeros((len(sub_group_list), 5))
        add_bias_sp = (len(sub_group_list), 5)

        # divide the group into many exposures
        for seq_tag, group_tag in enumerate(sub_group_list):
            # select individual group
            idx_group = group_label == group_tag
            sub_data = area_data[idx_group]
            group_ra = sub_data[:,5]
            group_dec = sub_data[:,6]

            sub_expos_labels = expos_labels[idx_group]

            group_ra_min, group_ra_max = group_ra.min(), group_ra.max()
            group_dec_min, group_dec_max = group_dec.min(), group_dec.max()

            group_ra_bin, group_ra_bin_num = hk_tool_box.set_min_bin(group_ra_min, group_ra_max, expos_field_width)
            group_dec_bin, group_dec_bin_num = hk_tool_box.set_min_bin(group_dec_min, group_dec_max, expos_field_width)

            # additive bias test
            # it should be very small
            gh1, gh1_sig = hk_FQlib.find_shear_cpp(sub_data[:, 0], sub_data[:, 2] + sub_data[:, 3], 20)[:2]
            gh2, gh2_sig = hk_FQlib.find_shear_cpp(sub_data[:, 1], sub_data[:, 2] - sub_data[:, 3], 20)[:2]
            # if the signal is significant than 1.5\sigma, it should be corrected
            if numpy.abs(gh1/gh1_sig) > 1.5:
                pass
                # sub_data[:, 0] = sub_data[:, 0] - gh1*(sub_data[:, 2] + sub_data[:, 3])
            if numpy.abs(gh2 / gh2_sig) > 1.5:
                pass
                # sub_data[:, 1] = sub_data[:, 1] - gh2*(sub_data[:, 2] - sub_data[:, 3])
            add_bias[seq_tag] = group_tag, gh1, gh1_sig,gh2, gh2_sig

            count = 0
            for i in range(group_ra_bin_num):
                idx1 = group_ra >= group_ra_bin[i]
                idx2 = group_ra < group_ra_bin[i+1]

                for j in range(group_dec_bin_num):
                    idx3 = group_dec >= group_dec_bin[j]
                    idx4 = group_dec < group_dec_bin[j+1]

                    idx_sub = idx1 & idx2 & idx3 & idx4

                    src_num = idx_sub.sum()

                    if src_num > 0:
                        expos_data = sub_data[idx_sub]
                        sub_redshift = expos_data[:, 8]
                        redshift_label = hk_tool_box.get_bin_label(sub_redshift, redshift_bin,redshift_bin_num)

                        expos_ra_center = (group_ra_bin[i] + group_ra_bin[i+1]) / 2
                        expos_dec_center = (group_dec_bin[j] + group_dec_bin[j+1]) / 2
                        cos_expos_dec_center = numpy.cos(expos_dec_center / deg2arcmin * deg2rad)

                        expos_pos = numpy.array([expos_ra_center, expos_dec_center, dra, ddec,
                                                numpy.sqrt((dra * cos_expos_dec_center) ** 2 + ddec ** 2), cos_expos_dec_center],
                                               dtype=numpy.float32)

                        expos_name = "%d-%d"%(group_tag+ncent_before, count)
                        expos_path = group_cata_path + "/%s.hdf5"%expos_name
                        h5f_expos = h5py.File(expos_path, "w")
                        h5f_expos["/pos"] = expos_pos
                        h5f_expos["/redshift_label"] = redshift_label
                        h5f_expos["/data"] = expos_data
                        h5f_expos["/expos_label"] = sub_expos_labels[idx_sub]
                        h5f_expos["/group_label"] = numpy.array([group_tag+ncent_before],dtype=numpy.intc)
                        h5f_expos.close()

                        expos_avail_sub.append("%s\t%s\t%d\t%f\t%f\t%f\t%f\t%f\t%f\n"
                                              % (expos_path, expos_name, src_num, expos_pos[0], expos_pos[1],
                                                 expos_pos[2], expos_pos[3], expos_pos[4], expos_pos[5]))
                        count += 1
                        expos_count += 1
        comm.Barrier()

        total_add_bias_sp = comm.gather(add_bias_sp, root=0)
        comm.Barrier()
        if rank > 0:
            # !!!! remember the data type, MPI.DOUBLE, MPI.FLOAT, ...
            # or it will raise an error, Keyerror
            comm.Send([add_bias, MPI.DOUBLE], dest=0, tag=rank)
        else:
            # receive the data from other CPUs
            # !!!! the start points is 1 in range() not 0
            for procs in range(1, cpus):
                # prepare a buffer for the data, the shape must be the same
                # with that of what the other CPUs send, you have collected them in 'data_sps'
                recvs = numpy.empty(total_add_bias_sp[procs], dtype=numpy.double)
                # receive it using the buffer,
                comm.Recv(recvs, source=procs, tag=procs)
                add_bias = numpy.row_stack((add_bias, recvs))

            h5f = h5py.File(result_cata_path + "/bias_test.hdf5", "r+")
            # print(results)
            h5f["/%d" % sub_area] = add_bias
            h5f.close()
        comm.Barrier()

    comm.Barrier()
    expos_count_list = comm.gather(expos_count, root=0)
    expos_avail_list = comm.gather(expos_avail_sub, root=0)

    if rank == 0:

        buffer_expo = []
        for fns in expos_avail_list:
            buffer_expo.extend(fns)

        with open(group_cata_path + "/source_list.dat", "w") as f:
            f.writelines(buffer_expo)

        log_inform = "%d exposures\n"%len(buffer_expo)
        with open("kmeans_log.dat", "w") as f:
            f.writelines(log_inform)
        print(log_inform)
        print(expos_count_list)
        print(sum(expos_count_list))
    comm.Barrier()
