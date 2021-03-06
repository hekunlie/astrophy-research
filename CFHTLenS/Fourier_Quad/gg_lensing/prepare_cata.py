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
bin_st, bin_ed = 0.04, 15
separation_bin = tool_box.set_bin_log(bin_st, bin_ed, sep_bin_num+1).astype(numpy.float32)

# bin number for ra & dec of each exposure
deg2arcmin = 60
deg2rad = numpy.pi/180

# position in CFHTLens catalog
ra_idx = 0
dec_idx = 1


redshift_idx = 10


# chi guess bin for PDF_SYM
chi_guess_num = 200
num_p = int(chi_guess_num/2)

delta_sigma_guess_bin_p = tool_box.set_bin_log(0.01, 500, num_p).astype(numpy.float64)
# delta_sigma_guess_bin_p = tool_box.set_bin_log(0.01, 130, num_p).astype(numpy.float64)
# delta_sigma_guess_bin_p = numpy.linspace(0.01, 150, num_p).astype(numpy.float64)

delta_sigma_guess = numpy.zeros((chi_guess_num, ), dtype=numpy.float64)
delta_sigma_guess[:num_p] = -delta_sigma_guess_bin_p
delta_sigma_guess[num_p:] = delta_sigma_guess_bin_p
delta_sigma_guess = numpy.sort(delta_sigma_guess)

gt_guess_num = 100
num_p = int(gt_guess_num/2)
# tan_shear_guess_bin_p = tool_box.set_bin_log(0.0005, 0.1, num_p).astype(numpy.float64)
tan_shear_guess_bin_p = numpy.linspace(0.0001, 0.1, num_p).astype(numpy.float64)

tan_shear_guess = numpy.zeros((gt_guess_num, ), dtype=numpy.float64)
tan_shear_guess[:num_p] = -tan_shear_guess_bin_p
tan_shear_guess[num_p:] = tan_shear_guess_bin_p
tan_shear_guess = numpy.sort(tan_shear_guess)



mg_bin_num = 20

hist2d_mg_num = 3000
hist2d_mg_num2 = int(hist2d_mg_num/2)

# star number on each chip
nstar_idx = 21
nstar_thresh = 12

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
Z_B_idx = 10
Z_B_MIN_idx = 38
Z_B_MAX_idx = 39
odd_idx = 40
odd_thresh = 0.5


#
# fourier_cata_path = "/coma/hklee/CFHT/CFHT_cat_Oct_11_2020"
# result_cata_path = "/coma/hklee/CFHT/correlation/cata"

fourier_cata_path = "/home/hklee/work/CFHT/CFHT_cat_4_20_2021"
result_cata_path = "/home/hklee/work/CFHT/gg_lensing/cata"
foreground_path_ori = "/home/hklee/work/catalog/cmass"

# fourier_cata_path = "/home/hklee/work/CFHT/CFHT_cat_Dec_17_2020_smoothed"
# result_cata_path = "/home/hklee/work/CFHT/gg_lensing/cata"

cmd = argv[1]



if cmd == "prepare_foreground":
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    cent_num = 200

    stack_file_path = foreground_path_ori + "/foreground.hdf5"

    files = ["cmass_w1.hdf5"]#, "cmass_w2.hdf5","cmass_w3.hdf5","cmass_w4.hdf5"]

    if rank == 0:
        if not os.path.exists(result_cata_path + "/foreground"):
            os.makedirs(result_cata_path + "/foreground")
        if not os.path.exists(result_cata_path + "/background"):
            os.makedirs(result_cata_path + "/background")

        src_num = numpy.zeros((len(files), ), dtype=numpy.intc)
        src_st = numpy.zeros((len(files), ), dtype=numpy.intc)

        for i, fn in enumerate(files):
            h5f = h5py.File(foreground_path_ori + "/" + fn, "r")
            data = h5f["/Z"][()]
            h5f.close()
            src_num[i] = data.shape[0]
            src_st[i] = src_num[:i].sum()

        total_num = src_num.sum()

        total_data = numpy.zeros((total_num, 6), dtype=numpy.float32)
        for i, fn in enumerate(files):
            st, ed = src_st[i], src_st[i] + src_num[i]
            h5f = h5py.File(foreground_path_ori + "/" + fn, "r")
            total_data[st:ed,0] = h5f["/RA"][()]
            total_data[st:ed,1] = h5f["/DEC"][()]
            total_data[st:ed,3] = h5f["/Z"][()]
            h5f.close()

        total_data[:, 2] = numpy.cos(total_data[:, 1]/180*numpy.pi)
        total_data[:, 4] = cosmos.comoving_distance(total_data[:, 3]).value*H0/100

        print(" %d galaxies" % (total_num))
        # Kmeans method for classification for jackknife
        t1 = time.time()

        total_data[:,5] = KMeans(n_clusters=cent_num, random_state=numpy.random.randint(1, 100000)).fit_predict(total_data[:,:2])

        h5f = h5py.File(stack_file_path, "w")
        h5f["/data"] = total_data
        h5f.close()

        t2 = time.time()
        print("Time: %.2f sec. %d galaxies"%(t2-t1, total_num))
    comm.Barrier()

    if rank > 0:
        h5f = h5py.File(stack_file_path, "r")
        total_data = h5f["/data"][()]
        h5f.close()

    # assign the source into the artificial exposures
    min_src_num = 100

    expos_avail_sub = []
    expos_count = 0

    group_list = [i for i in range(cent_num)]

    sub_group_list = tool_box.alloc(group_list, cpus)[rank]
    group_label = total_data[:,5].astype(dtype=numpy.intc)

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
            nums_distrib = tool_box.alloc([1 for i in range(ground_src_num)], m)
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


elif cmd == "prepare_background":

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

            ra = src_data[:, ra_idx]
            dec = src_data[:, dec_idx]

            # find the center and the 4 corners
            ra_min, ra_max = ra.min(), ra.max()
            ra_center = (ra_min + ra_max) / 2
            dra = (ra_max - ra_min) / 2
            dec_min, dec_max = dec.min(), dec.max()
            dec_center = (dec_min + dec_max) / 2
            cos_dec_center = numpy.cos(dec_center * deg2rad)
            ddec = (dec_max - dec_min) / 2

            expo_pos = numpy.array([ra_center, dec_center, cos_dec_center,
                                     numpy.sqrt((dra * cos_dec_center) ** 2 + ddec ** 2)], dtype=numpy.float32)

            # G1, G2, N, U, V, RA, DEC, COS(DEC), Z, Z_ERR, COMOVING DISTANCE
            dst_data = numpy.zeros((src_num, 11), dtype=numpy.float32)

            dst_data[:, 0] = src_data[:, mg1_idx]
            dst_data[:, 1] = src_data[:, mg2_idx]
            dst_data[:, 2] = src_data[:, mn_idx]
            # the signs of U & V are different with that in the paper
            dst_data[:, 3] = -src_data[:, mu_idx]
            dst_data[:, 4] = -src_data[:, mv_idx]

            dst_data[:, 5] = src_data[:, ra_idx]
            dst_data[:, 6] = src_data[:, dec_idx]
            dst_data[:, 7] = numpy.cos(dst_data[:, 6] * deg2rad)
            dst_data[:, 8] = src_data[:, Z_B_idx]
            dst_data[:, 9] = (src_data[:, Z_B_MAX_idx] - src_data[:,Z_B_MIN_idx])/2
            dst_data[:, 10] = cosmos.comoving_distance(dst_data[:,8]).value*H0/100 # in unit of Mpc/h

            expo_dst_path = result_cata_path + "/background/%s" %expo_name
            h5f_dst = h5py.File(expo_dst_path,"w")
            h5f_dst["/pos"] = expo_pos
            h5f_dst["/data"] = dst_data
            h5f_dst.close()

            expo_avail_sub.append("%s\t%s\t%d\t%f\t%f\t%f\t%f\n"
                                   %(expo_dst_path, expo_name.split(".")[0], src_num, expo_pos[0],expo_pos[1], expo_pos[2],
                                     expo_pos[3]))
        # exit()
    comm.Barrier()

    expo_avail_list = comm.gather(expo_avail_sub, root=0)
    exception_collection = comm.gather(exception_sub, root=0)

    if rank == 0:
        buffer_expo = []
        for fns in expo_avail_list:
            buffer_expo.extend(fns)

        with open(result_cata_path + "/background/background_source_list.dat", "w") as f:
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


if cmd == "prepare_pdf":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

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

        # G bins for tangential shear calculation
        mg_bin = tool_box.set_bin(src_data[:, 0], mg_bin_num, 100000)

        # G_t bins for \Delta\Sigma(R) calculation
        Gts_total_num = 1000000
        Gts = numpy.zeros((Gts_total_num,), dtype=numpy.float32)
        NU = numpy.zeros((Gts_total_num,), dtype=numpy.float32)
        Gt_num = 0

        with open(result_cata_path + "/background/background_source_list.dat", "r") as f:
            back_contents = f.readlines()
        back_expo_num = len(back_contents)
        back_expo_cent = numpy.zeros((back_expo_num, 2),dtype=numpy.float32)
        back_expo_path = []
        back_expo_labels = numpy.arange(0,back_expo_num)
        for back_c in range(back_expo_num):
            cc = back_contents[back_c].split("\n")[0].split("\t")
            back_expo_path.append(cc[0])
            back_expo_cent[back_c] = float(cc[3]),float(cc[4])
        back_expo_skypos = SkyCoord(ra=back_expo_cent[:,0]*units.deg, dec=back_expo_cent[:,1]*units.deg,frame="fk5")


        with open(result_cata_path + "/foreground/foreground_source_list.dat", "r") as f:
            fore_contents = f.readlines()

        foreground_num = len(fore_contents)
        for fore_c in fore_contents:
            cc = fore_c.split("\n")[0].split("\t")
            fore_expo_path = cc[0]
            h5f = h5py.File(fore_expo_path, "r")
            fore_data = h5f["/data"][()]
            h5f.close()

            fore_c_num = fore_data.shape[0]

            for ic in range(fore_c_num):
                if Gt_num >= Gts_total_num:
                    break
                ic_z = fore_data[ic,3]
                ic_com_dist = fore_data[ic, 4]
                ic_skypos = SkyCoord(ra=fore_data[ic,0]*units.deg, dec=fore_data[ic,1]*units.deg, frame="fk5")
                separation_angle = ic_skypos.separation(back_expo_skypos).deg
                idx = separation_angle < 2
                if idx.sum() > 0:
                    for ib in back_expo_labels[idx]:

                        if Gt_num >= Gts_total_num:
                            break

                        h5f = h5py.File(back_expo_path[ib], "r")
                        back_data = h5f["/data"][()]
                        h5f.close()

                        idxz = back_data[:,8] > ic_z + 0.1
                        sub_num = idxz.sum()
                        if sub_num > 0:
                            sub_data = back_data[idxz]
                            ib_skypos = SkyCoord(ra=sub_data[:,5]*units.deg, dec=sub_data[:,6]*units.deg,frame="fk5")
                            position_ang = 2*ic_skypos.position_angle(ib_skypos).radian

                            mgt = sub_data[:,0]*numpy.cos(position_ang) - sub_data[:,1]*numpy.sin(position_ang)
                            mnu = sub_data[:, 2] + sub_data[:, 3]*numpy.cos(2*position_ang) - sub_data[:,4]*numpy.sin(2 * position_ang)

                            # 1662895.2007121066 = c^2/4/pi/G  [M_sum/pc] /10^6 [h/pc]
                            mgt = mgt*sub_data[:,10]/ic_com_dist/(sub_data[:,10]-ic_com_dist)*1662895.2007121066
                            if Gt_num + sub_num > Gts_total_num:
                                Gts[Gt_num: Gts_total_num] = mgt[Gts_total_num-Gt_num]
                                NU[Gt_num: Gts_total_num] = mnu[Gts_total_num - Gt_num]
                                Gt_num = Gts_total_num
                                break
                            else:
                                Gts[Gt_num: Gt_num + sub_num] = mgt
                                NU[Gt_num: Gt_num + sub_num] = mnu
                                Gt_num += sub_num

        mg_sigma_bin = tool_box.set_bin(Gts, mg_bin_num, 1000000)
        print(tan_shear_guess)
        print(delta_sigma_guess)
        print(mg_bin)
        print(mg_sigma_bin)

        hist2d_mg_bin = numpy.zeros((hist2d_mg_num + 1,))
        hist2d_mnu_bin = numpy.zeros((hist2d_mg_num + 1,))

        hist2d_mg_bin[:hist2d_mg_num2] = -tool_box.set_bin_log(0.01, Gts.max() * 2000, hist2d_mg_num2)
        hist2d_mg_bin[hist2d_mg_num2 + 1:] = tool_box.set_bin_log(0.01, Gts.max() * 2000, hist2d_mg_num2)

        hist2d_mnu_bin[:hist2d_mg_num2] = -tool_box.set_bin_log(0.01, NU.max() * 2000, hist2d_mg_num2)
        hist2d_mnu_bin[hist2d_mg_num2 + 1:] = tool_box.set_bin_log(0.01, NU.max() * 2000, hist2d_mg_num2)

        hist2d_mg_bin = numpy.sort(hist2d_mg_bin)
        hist2d_mnu_bin = numpy.sort(hist2d_mnu_bin)


        h5f = h5py.File(result_cata_path + "/pdf_inform.hdf5", "w")

        h5f["/hist2d_mg_sigma_bin"] = hist2d_mg_bin.astype(dtype=numpy.float32)
        h5f["/hist2d_mn_sigma_bin"] = hist2d_mnu_bin.astype(dtype=numpy.float32)

        h5f["/mg_sigma_bin"] = numpy.array(mg_sigma_bin, dtype=numpy.float32)
        h5f["/mg_gt_bin"] = numpy.array(mg_bin, dtype=numpy.float32)
        h5f["/gt_guess"] = tan_shear_guess
        h5f["/delta_sigma_guess"] = delta_sigma_guess
        h5f["/separation_bin"] = separation_bin
        h5f["/cosmological_params"] = numpy.array([H0, omega_m0], dtype=numpy.float32)

        h5f.close()

    comm.Barrier()
