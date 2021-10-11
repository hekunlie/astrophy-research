import os
my_home = os.popen("echo $HK_MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import shutil
import h5py
import numpy
from mpi4py import MPI
import hk_tool_box
import warnings
from sklearn.cluster import KMeans
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
import time


warnings.filterwarnings('error')


# parameters
# cosmological parameters
omega_m0 = 0.31
H0 = 67.5
cosmos = FlatLambdaCDM(H0, omega_m0)

# separation bin, comoving or angular diameter distance in unit of Mpc/h
sep_bin_num = 13
bin_st, bin_ed = 0.05, 20
separation_bin = hk_tool_box.set_bin_log(bin_st, bin_ed, sep_bin_num+1).astype(numpy.float32)

# bin number for ra & dec of each exposure
deg2arcmin = 60
deg2rad = numpy.pi/180


# chi guess bin for PDF_SYM
delta_sigma_guess_num = 100
num_m = 50
num_p = delta_sigma_guess_num - num_m

delta_sigma_guess = numpy.zeros((delta_sigma_guess_num, ), dtype=numpy.float64)

delta_sigma_guess_bin_p = hk_tool_box.set_bin_log(0.001, 400, num_p).astype(numpy.float64)

delta_sigma_guess[:num_m] = -hk_tool_box.set_bin_log(0.001, 400, num_m).astype(numpy.float64)
delta_sigma_guess[num_m:] = hk_tool_box.set_bin_log(0.001, 400, num_p).astype(numpy.float64)
delta_sigma_guess = numpy.sort(delta_sigma_guess)


gt_guess_num = 400
num_p = int(gt_guess_num/2)


# tan_shear_guess_bin_p = tool_box.set_bin_log(0.0005, 0.2, num_p).astype(numpy.float64)
tan_shear_guess_bin_p = numpy.linspace(0.0001, 0.1, num_p).astype(numpy.float64)


tan_shear_guess = numpy.zeros((gt_guess_num, ), dtype=numpy.float64)
tan_shear_guess[:num_p] = -tan_shear_guess_bin_p
tan_shear_guess[num_p:] = tan_shear_guess_bin_p
tan_shear_guess = numpy.sort(tan_shear_guess)


mg_bin_num = 10

hist2d_mg_num = 4
hist2d_mg_num2 = int(hist2d_mg_num/2)

# position in DECALS catalog
ra_idx = 0
dec_idx = 1

# star number on each chip
nstar_idx = 6
nstar_thresh = 20

# selection function
flux2_alt_idx = 7
flux2_alt_thresh = 3


# field distortion
gf1_idx = 8
gf2_idx = 9
gf1_thresh = 0.0015
gf2_thresh = 0.0015

# shear estimators
mg1_idx = 10
mg2_idx = 11
mn_idx = 12
mu_idx = 13
mv_idx = 14

# # delta ra dec
# # the pixel scale is 0.262 arcsec/pix
# # 1 pix corresponds to 7.28 * 10^{-5}
# dist_idx = 13
# dist_thresh = 10**(-4)
#
# # PhotoZ
Zp_idx = 4 # photo Z
# Zs_idx = 17 # spectral Z


# fourier_cata_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/DECALS"
# result_cata_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/DECALS/gg_lensing"
# foreground_path_ori = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/Yang_group"

fourier_cata_path = "/home/hklee/work/DECALS/DECALS_shear_catalog_v210729"
result_cata_path = "/home/hklee/work/DECALS/DECALS_shear_catalog_v210729/gg_lensing/cata"
# foreground_path_ori = "/home/hklee/work/catalog/Yang_group"
foreground_path_ori = "/home/hklee/work/catalog/SDSS"
# foreground_path_ori = "/home/hklee/work/catalog/Jesse_cata/hdf5"
fourier_avail_expo_path = fourier_cata_path + "/cat_inform/exposure_avail_rz_band.dat"


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

####### foreground selection #######
fore_richness_idx = 0
fore_ra_idx = 1
fore_dec_idx = 2
fore_z_idx = 3
fore_mass_idx = 4  # log M

fore_richness_thresh = 4
fore_z_min = float(argv[2])#0.3
fore_z_max = float(argv[3])#0.4
fore_mass_min = 1#float(argv[4])#13.5
fore_mass_max = 20#float(argv[5])#13
if rank == 0:
    log_inform = "Foreground selection: richness>=%d, " \
                 "Z: %.2f~%.2f, Mass: %.2f~%.2f"%(fore_richness_thresh, fore_z_min,fore_z_max,fore_mass_min,fore_mass_max)
    print(log_inform)


stack_file_path = foreground_path_ori + "/foreground.hdf5"
rand_stack_file_path = foreground_path_ori + "/foreground_rand.hdf5"

# files = ["DESI_NGC_group_DECALS_overlap.hdf5","DESI_SGC_group_DECALS_overlap.hdf5"]
files = ["lowz_DECALS_overlap.hdf5","cmass_DECALS_overlap.hdf5"]
# files = ["lowz_DECALS_overlap.hdf5"]
# files = [argv[2]]

# for kmeans to build jackknife labels
cent_num = 200

if rank == 0:
    if os.path.exists(result_cata_path + "/foreground"):
        shutil.rmtree(result_cata_path + "/foreground")
    os.makedirs(result_cata_path + "/foreground")

    for i, fn in enumerate(files):
        fore_path = foreground_path_ori + "/" + fn
        print("\nRead %s\n"%fore_path)

        h5f = h5py.File(fore_path, "r")
        temp = h5f["/data"][()]
        h5f.close()

        if i == 0:
            data_src = temp
        else:
            data_src = numpy.row_stack((data_src, temp))

    # foreground selection
    idx_z1 = data_src[:,fore_z_idx] >= fore_z_min
    idx_z2 = data_src[:,fore_z_idx] < fore_z_max
    # idx_r = data_src[:, fore_richness_idx] >= fore_richness_thresh
    # idx_m1 = data_src[:,fore_mass_idx] >= fore_mass_min
    # idx_m2 = data_src[:,fore_mass_idx] < fore_mass_max
    idx = idx_z1 & idx_z2 #& idx_r & idx_m1 & idx_m2

    total_num = idx.sum()
    # total_num = data_src.shape[0]

    total_data = numpy.zeros((total_num, 5), dtype=numpy.float32)
    for i, fn in enumerate(files):
        h5f = h5py.File(foreground_path_ori + "/" + fn, "r")
        total_data[:,0] = data_src[:,fore_ra_idx][idx]
        total_data[:,1] = data_src[:,fore_dec_idx][idx]
        total_data[:,3] = data_src[:,fore_z_idx][idx]
        h5f.close()

    total_data[:, 2] = numpy.cos(total_data[:, 1]*deg2rad)
    total_data[:, 4] = cosmos.comoving_distance(total_data[:, 3]).value*H0/100

    # the random sample
    rand_total_data = numpy.zeros((total_num, 5), dtype=numpy.float32)
    h5f = h5py.File(foreground_path_ori + "/rand_fore_candidate.hdf5", "r")
    temp = h5f["/data"][()]
    h5f.close()
    rand_total_data[:, 0] = temp[:total_num, fore_ra_idx]
    rand_total_data[:, 1] = temp[:total_num, fore_dec_idx]
    rand_total_data[:, 2] = numpy.cos(rand_total_data[:, 1]*deg2rad)
    # use the len z for the random points
    rand_total_data[:, 3] = total_data[:, 3]
    rand_total_data[:, 4] = total_data[:, 4]


    print(" %d galaxies" % (total_num))
    # Kmeans method for classification for jackknife
    t1 = time.time()

    rs = numpy.random.randint(1, 100000)

    group_label = KMeans(n_clusters=cent_num, random_state=rs).fit_predict(total_data[:,:2])

    # random sample
    rand_group_label = KMeans(n_clusters=cent_num, random_state=rs).fit_predict(rand_total_data[:,:2])

    idx_1 = group_label >= 0
    idx_2 = group_label < cent_num
    idx_label = idx_1 & idx_2
    if idx_label.sum() != total_num:
        print("KMeans wrong!!! ",idx_label.sum(),total_num)

    t2 = time.time()

    h5f = h5py.File(stack_file_path, "w")
    h5f["/data"] = total_data
    h5f["/group_label"] = group_label.astype(dtype=numpy.intc)
    h5f.close()

    # random sample
    h5f = h5py.File(rand_stack_file_path, "w")
    h5f["/data"] = rand_total_data
    h5f["/group_label"] = rand_group_label.astype(dtype=numpy.intc)
    h5f.close()

    print("Time: %.2f sec. %d galaxies"%(t2-t1, total_num))
comm.Barrier()

if rank > 0:
    h5f = h5py.File(stack_file_path, "r")
    total_data = h5f["/data"][()]
    group_label = h5f["/group_label"][()]
    h5f.close()

    # random sample
    h5f = h5py.File(rand_stack_file_path, "r")
    rand_total_data = h5f["/data"][()]
    rand_group_label = h5f["/group_label"][()]
    h5f.close()

comm.Barrier()


# assign the source into the artificial exposures
min_src_num = 50

expos_avail_sub = []
expos_count = 0

group_list = [i for i in range(cent_num)]

sub_group_list = hk_tool_box.alloc(group_list, cpus)[rank]

# divide the group into many exposures
for group_tag in sub_group_list:

    # select individual group
    idx_group = group_label == group_tag
    sub_data = total_data[idx_group]

    ground_src_num = idx_group.sum()

    if ground_src_num <= min_src_num:
        expos_name = "%d-0" %group_tag
        expos_path = result_cata_path + "/foreground/%s.hdf5" % expos_name
        h5f_expos = h5py.File(expos_path, "w")
        h5f_expos["/data"] = sub_data
        h5f_expos.close()

        expos_avail_sub.append("%s\t%s\t%d\t%d\n"
                               % (expos_path, expos_name, ground_src_num, group_tag))
        expos_count += 1
    else:
        m, n = divmod(ground_src_num, min_src_num)
        nums_distrib = hk_tool_box.alloc([1 for i in range(ground_src_num)], m)
        nums = [sum(nums_distrib[i]) for i in range(m)]
        nums_st = [sum(nums[:i]) for i in range(m)]
        for count in range(m):
            expos_name = "%d-%d" % (group_tag, count)
            expos_path = result_cata_path + "/foreground/%s.hdf5" % expos_name
            h5f_expos = h5py.File(expos_path, "w")
            h5f_expos["/data"] = sub_data[nums_st[count]: nums_st[count]+nums[count]]
            h5f_expos.close()

            expos_avail_sub.append("%s\t%s\t%d\t%d\n"
                                   % (expos_path, expos_name, nums[count], group_tag))
            expos_count += 1

comm.Barrier()
expos_count_list = comm.gather(expos_count, root=0)
expos_avail_list = comm.gather(expos_avail_sub, root=0)

if rank == 0:

    buffer_expo = []
    for fns in expos_avail_list:
        buffer_expo.extend(fns)

    with open(result_cata_path + "/foreground/foreground_source_list.dat", "w") as f:
        f.writelines(buffer_expo)

    log_inform = "%d exposures\n"%len(buffer_expo)
    with open("kmeans_log.dat", "w") as f:
        f.writelines(log_inform)
    print(log_inform)
    print(expos_count_list)
    print(sum(expos_count_list))
comm.Barrier()

# the random sample

rand_expos_avail_sub = []
rand_expos_count = 0

rand_group_list = [i for i in range(cent_num)]

sub_rand_group_list = hk_tool_box.alloc(rand_group_list, cpus)[rank]

# divide the group into many exposures
for group_tag in sub_group_list:

    # select individual group
    idx_group = group_label == group_tag
    sub_data = rand_total_data[idx_group]

    ground_src_num = idx_group.sum()

    if ground_src_num <= min_src_num:
        expos_name = "%d-0" %group_tag
        expos_path = result_cata_path + "/foreground/%s.hdf5" % expos_name
        h5f_expos = h5py.File(expos_path, "w")
        h5f_expos["/data"] = sub_data
        h5f_expos.close()

        expos_avail_sub.append("%s\t%s\t%d\t%d\n"
                               % (expos_path, expos_name, ground_src_num, group_tag))
        expos_count += 1
    else:
        m, n = divmod(ground_src_num, min_src_num)
        nums_distrib = hk_tool_box.alloc([1 for i in range(ground_src_num)], m)
        nums = [sum(nums_distrib[i]) for i in range(m)]
        nums_st = [sum(nums[:i]) for i in range(m)]
        for count in range(m):
            expos_name = "%d-%d" % (group_tag, count)
            expos_path = result_cata_path + "/foreground/%s.hdf5" % expos_name
            h5f_expos = h5py.File(expos_path, "w")
            h5f_expos["/data"] = sub_data[nums_st[count]: nums_st[count]+nums[count]]
            h5f_expos.close()

            expos_avail_sub.append("%s\t%s\t%d\t%d\n"
                                   % (expos_path, expos_name, nums[count], group_tag))
            expos_count += 1

comm.Barrier()
expos_count_list = comm.gather(expos_count, root=0)
expos_avail_list = comm.gather(expos_avail_sub, root=0)

if rank == 0:

    buffer_expo = []
    for fns in expos_avail_list:
        buffer_expo.extend(fns)

    with open(result_cata_path + "/foreground/foreground_source_list.dat", "w") as f:
        f.writelines(buffer_expo)

    log_inform = "%d exposures\n"%len(buffer_expo)
    with open("kmeans_log.dat", "w") as f:
        f.writelines(log_inform)
    print(log_inform)
    print(expos_count_list)
    print(sum(expos_count_list))
comm.Barrier()





