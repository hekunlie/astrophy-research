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
from astropy.cosmology import FlatLambdaCDM
import time
import prepare_tools

warnings.filterwarnings('error')


# parameters
# cosmological parameters
omega_m0 = 0.27
H0 = 67.5
cosmos = FlatLambdaCDM(H0, omega_m0)

area_num = 4
# separation bin
sep_bin_num = 5
separation_bin = tool_box.set_bin_log(0.8, 60, sep_bin_num+1).astype(numpy.float32)

# bin number for ra & dec of each exposure
deg2arcmin = 60
deg2rad = numpy.pi/180

grid_size = 15 #arcmin

# position in CFHTLens catalog
ra_idx = 0
dec_idx = 1


redshift_idx = 10


# chi guess bin for PDF_SYM
chi_guess_num = 150
num_p = int(6/8*chi_guess_num)
num_m = chi_guess_num - num_p

delta_sigma_guess_bin_p = numpy.linspace(0.1, 600, num_p).astype(numpy.float32)
delta_sigma_guess_bin_m = numpy.linspace(0.1, 200, num_m).astype(numpy.float32)

delta_sigma_guess = numpy.zeros((chi_guess_num, ), dtype=numpy.float32)
delta_sigma_guess[:num_m] = -delta_sigma_guess_bin_m
delta_sigma_guess[num_m:] = delta_sigma_guess_bin_p
delta_sigma_guess = numpy.sort(delta_sigma_guess)

tan_shear_guess = numpy.linspace(-0.05, 0.1, chi_guess_num, dtype=numpy.float32)


mg_bin_num = 10

# star number on each chip
nstar_idx = 21
nstar_thresh = 20

# selection function
flux2_alt_idx = 28
flux2_alt_thresh = 2.5

# selection on the bad pixels
imax_idx = 22
jmax_idx = 23
imax_thresh = 48
jmax_thresh = 48

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
odd_idx = 38
odd_thresh = 0.5
redshift_e_idx = 39
pz_sum_idx = 40
# pz_sum_thresh = 0.001

#
# fourier_cata_path = "/coma/hklee/CFHT/CFHT_cat_Oct_11_2020"
# result_cata_path = "/coma/hklee/CFHT/correlation/cata"

# fourier_cata_path = "/mnt/perc/hklee/CFHT/CFHT_cat_Dec_17_2020_smoothed"
# result_cata_path = "/mnt/perc/hklee/CFHT/correlation/cata"

fourier_cata_path = "/home/hklee/work/CFHT/CFHT_cat_Dec_17_2020_smoothed"
result_cata_path = "/home/hklee/work/CFHT/gg_lensing/cata"

cmd = argv[1]

if cmd == "prepare_pdf":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()
    # rank = 0
    # cpus = 1

    # for correlation calculation
    if rank == 0:

        field_name = []
        with open(fourier_cata_path + "/cat_inform/exposure_avail.dat", "r") as f:
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

        mg_bin = tool_box.set_bin(src_data[:, 0], mg_bin_num, 100000)

        h5f = h5py.File(result_cata_path + "/pdf_inform.hdf5", "w")

        h5f["/mg_bin"] = mg_bin
        h5f["/g_guess"] = tan_shear_guess
        h5f["/delta_sigma_guess"] = delta_sigma_guess
        h5f["/separation_bin"] = separation_bin
        h5f["/cosmological_params"] = numpy.array([H0, omega_m0], dtype=numpy.float32)
        h5f.close()

    comm.Barrier()


if cmd == "prepare_background":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    total_expos = []
    with open(fourier_cata_path + "/cat_inform/exposure_avail.dat", "r") as f:
        f_lines = f.readlines()
    for ff in f_lines:
        total_expos.append(ff.split("\n")[0])

    my_expos = tool_box.alloc(total_expos,cpus, method="order")[rank]
    expo_avail_sub = []
    exception_sub = []

    for cat_path in my_expos:
        expo_name = cat_path.split("/")[-1]

        h5f_src = h5py.File(cat_path,"r")
        data = h5f_src["/data"][()]
        h5f_src.close()

        # selection
        idx1 = data[:, nstar_idx] >= nstar_thresh
        idx2 = data[:, flux2_alt_idx] >= flux2_alt_thresh
        idx3 = numpy.abs(data[:, gf1_idx]) <= gf1_thresh
        idx4 = numpy.abs(data[:, gf2_idx]) <= gf2_thresh
        idx5 = data[:, imax_idx] < imax_thresh
        idx6 = data[:, jmax_idx] < jmax_thresh
        idx7 = data[:,odd_idx] >= odd_thresh
        idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx7

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
            # G1, G2, N, U, V, RA, DEC, COS(DEC), REDSHIFT, COM_DISTANCE
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
            dst_data[:, 9] = cosmos.comoving_distance(dst_data[:,8]).value*H0/100 # in unit of h/Mpc

            expo_dst_path = result_cata_path + "/%s" %expo_name
            h5f_dst = h5py.File(expo_dst_path,"w")
            h5f_dst["/pos"] = expo_pos
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

        with open(result_cata_path + "/background_source_list.dat", "w") as f:
            f.writelines(buffer_expo)

        exception_all = []
        for ec in exception_collection:
            exception_all.extend(ec)
        log_inform = "%d (%d) exposures\n"%(len(buffer_expo), len(total_expos))
        exception_all.append(log_inform)
        with open("log.dat", "w") as f:
            f.writelines(exception_all)
        print(log_inform)
    comm.Barrier()

elif cmd == "kmeans":
    # Kmeans method for classification for jackknife
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    t1 = time.time()

    # ncent = int(argv[2])
    # njobs = int(argv[3])
    # 200 sub-samples
    # ncent = [75, 34, 56, 35][rank]
    ra_bin = [1000, 5000, 10000, 15000, 30000]
    total_cent = 200

    if rank < len(ra_bin) - 1:

        h5f = h5py.File(result_cata_path + "/stack_data.hdf5", "r")
        data = h5f["/data"][()]
        expos_labels = h5f["/expos_label"][()]
        h5f.close()

        total_src_num = data.shape[0]

        ra_dec = data[:,5:7]
        src_num_each_area = numpy.zeros((area_num,),dtype=numpy.intc)
        for i in range(area_num):
            idx1 = ra_dec[:,0] >= ra_bin[i]
            idx2 = ra_dec[:,0] < ra_bin[i+1]
            idx = idx1 & idx2
            src_num_each_area[i] = idx.sum()

        ncents, ratio = prepare_tools.even_area(src_num_each_area, area_num, total_cent)
        if rank == 0:
            print(ncents)
            print(ratio)
        ncent = ncents[rank]
        idx1 = ra_dec[:, 0] >= ra_bin[rank]
        idx2 = ra_dec[:, 0] < ra_bin[rank + 1]
        idx = idx1 & idx2

        group_pred = KMeans(n_clusters=ncent, random_state=numpy.random.randint(1, 100000)).fit_predict(ra_dec[idx])

        h5f = h5py.File(result_cata_path + "/group_predict_%d.hdf5"%rank, "w")
        h5f["/data"] = group_pred
        h5f["/ra_dec"] = ra_dec[idx]
        h5f["/expos_labels"] = expos_labels[idx]
        h5f.close()

        t2 = time.time()
        print("%d Total sub-sample: %d (%d ~ %d). Time: %.2f sec."%(rank, ncent, group_pred.min(), group_pred.max(), t2-t1))
    comm.Barrier()

elif cmd == "segment":
    # assign the source into the artificial exposures

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    group_cata_path = result_cata_path + "/kmeans"
    if rank == 0:
        if not os.path.exists(group_cata_path):
            os.makedirs(group_cata_path)
    comm.Barrier()

    # area num
    are_num = 4
    # the width of each exposure, arcmin
    expos_field_width = 40
    dra = expos_field_width/2
    ddec = expos_field_width/2
    ra_bin = [1000, 5000, 10000, 15000, 30000]
    ncent = []

    expos_avail_sub = []
    expos_count = 0
    for sub_area in range(area_num):
        # read the group labels from Kmeans
        h5f = h5py.File(result_cata_path + "/group_predict_%d.hdf5"%sub_area,"r")
        group_label = h5f["/data"][()]
        expos_labels = h5f["/expos_labels"][()]
        h5f.close()

        # the catalog
        h5f = h5py.File(result_cata_path + "/stack_data.hdf5", "r")
        data = h5f["/data"][()]
        h5f.close()

        ra_dec = data[:, 5:7]
        idx1 = ra_dec[:, 0] >= ra_bin[sub_area]
        idx2 = ra_dec[:, 0] < ra_bin[sub_area + 1]
        idx_area = idx1 & idx2
        area_data = data[idx_area]

        group_num = group_label.max() + 1
        ncent_before = sum(ncent)
        ncent.append(group_num)

        group_list = [i for i in range(group_num)]

        sub_group_list = tool_box.alloc(group_list, cpus)[rank]

        # divide the group into many exposures
        for group_tag in sub_group_list:
            # select individual group
            idx_group = group_label == group_tag
            sub_data = area_data[idx_group]
            group_ra = sub_data[:,5]
            group_dec = sub_data[:,6]
            sub_expos_labels = expos_labels[idx_group]

            group_ra_min, group_ra_max = group_ra.min(), group_ra.max()
            group_dec_min, group_dec_max = group_dec.min(), group_dec.max()

            group_ra_bin, group_ra_bin_num = prepare_tools.set_min_bin(group_ra_min, group_ra_max, expos_field_width)
            group_dec_bin, group_dec_bin_num = prepare_tools.set_min_bin(group_dec_min, group_dec_max, expos_field_width)

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
                        redshift_label = prepare_tools.get_bin_label(sub_redshift, redshift_bin,redshift_bin_num)

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