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
import hk_healpy_tool


warnings.filterwarnings('error')


# parameters
# cosmological parameters
omega_m0 = 0.315
H0 = 67.5
cosmos = FlatLambdaCDM(H0, omega_m0)

# separation bin, comoving or angular diameter distance in unit of Mpc/h
sep_bin_num = 25
bin_st, bin_ed = 0.02, 100
separation_bin = hk_tool_box.set_bin_log(bin_st, bin_ed, sep_bin_num+1).astype(numpy.float32)

# bin number for ra & dec of each exposure
deg2arcmin = 60
deg2rad = numpy.pi/180


# chi guess bin for PDF_SYM
delta_sigma_guess_num = 50
num_m = 25
num_p = delta_sigma_guess_num - num_m

delta_sigma_guess = numpy.zeros((delta_sigma_guess_num, ), dtype=numpy.float64)

delta_sigma_guess_bin_p = hk_tool_box.set_bin_log(0.01, 400, num_p).astype(numpy.float64)
# delta_sigma_guess_bin_p = numpy.array([1000/1.3**i for i in range(num_p)])#
delta_sigma_guess[:num_m] = -delta_sigma_guess_bin_p
delta_sigma_guess[num_m:] = delta_sigma_guess_bin_p
delta_sigma_guess = numpy.sort(delta_sigma_guess)


gt_guess_num = 100
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
flux2_alt_thresh = 2


# field distortion
gf1_idx = 8
gf2_idx = 9
gf1_thresh = 0.0015
gf2_thresh = 0.0015
gf_thresh = 0.05

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

fourier_cata_path = "/home/hklee/work/DECALS/DECALS_v210729"
result_cata_path = "/home/hklee/work/DECALS/DECALS_v210729/gg_lensing/cata"
# foreground_path_ori = "/home/hklee/work/catalog/Yang_group"
# foreground_path_ori = "/home/hklee/work/catalog/SDSS"
# foreground_path_ori = "/home/hklee/work/catalog/Jesse_cata/hdf5"
# foreground_path_ori = "/home/hklee/work/catalog/Yang_group/DESI_CLUSTER_DR9"
foreground_path_ori = "/home/hklee/work/catalog/Haojie_cata"
fourier_avail_expo_path = fourier_cata_path + "/cat_inform/exposure_avail_rz_band.dat"

parent_path_pi2 = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/DECALS_v210729/gg_lensing/cata"


cmd = argv[1]

if cmd == "prepare_foreground":

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
    fore_z_min = 0#float(argv[2])#0.3
    fore_z_max = 10#float(argv[3])#0.4
    fore_mass_min = 12#float(argv[4])#13.5
    fore_mass_max = 20#float(argv[5])#13
    if rank == 0:
        log_inform = "H0: %.4f. Omega_m0: %.4f"%(H0, omega_m0)
        print(log_inform)

        log_inform = "Foreground selection: richness>=%d, " \
                     "Z: %.2f~%.2f, Mass: %.2f~%.2f"%(fore_richness_thresh, fore_z_min,fore_z_max,fore_mass_min,fore_mass_max)
        print(log_inform)


    stack_file_path = foreground_path_ori + "/foreground.hdf5"

    # files = ["DESI_NGC_group_DECALS_overlap.hdf5","DESI_SGC_group_DECALS_overlap.hdf5"]
    # files = ["lowz_DECALS_overlap.hdf5","cmass_DECALS_overlap.hdf5"]
    # files = ["lowz_DECALS_overlap.hdf5"]
    files = ["gal_jkf_dr9.hdf5"]

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
            group_label = h5f["/jkf_label_200"][()]
            # pix_ra_dec, pixel_count, src_pix_label = h5f["/pix_ra_dec"][()], h5f["/pix_count"][()], h5f["/data_pix_label"][()]
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
        idx = idx_z1 & idx_z2 #& idx_m1 & idx_m2#& idx_r & idx_m1 & idx_m2

        total_num = idx.sum()
        # total_num = data_src.shape[0]
        # total_num = 450000
        total_data = numpy.zeros((total_num, 8), dtype=numpy.float32)

        for i, fn in enumerate(files):
            h5f = h5py.File(foreground_path_ori + "/" + fn, "r")
            total_data[:, 0] = data_src[:,fore_ra_idx][idx][:total_num]
            total_data[:, 1] = data_src[:,fore_ra_idx][idx][:total_num]*deg2rad
            total_data[:, 2] = data_src[:,fore_dec_idx][idx][:total_num]
            total_data[:, 3] = data_src[:,fore_dec_idx][idx][:total_num]*deg2rad
            total_data[:, 6] = data_src[:,fore_z_idx][idx][:total_num]
            h5f.close()

        total_data[:, 4] = numpy.cos(total_data[:, 3])
        total_data[:, 5] = numpy.sin(total_data[:, 3])
        total_data[:, 7] = cosmos.comoving_distance(total_data[:, 6]).value*H0/100

        print(" %d galaxies" % (total_num))
        # Kmeans method for classification for jackknife
        t1 = time.time()

        rs = numpy.random.randint(1, 100000)

        # eff_pix_num = len(pixel_count)
        # group_label, group_label_pixel = hk_healpy_tool.kmeans_pix(pix_ra_dec, pixel_count, src_pix_label, eff_pix_num,
        #                                                            cent_num, rs)

        # if not os.path.exists(foreground_path_ori + "/jkf_label.hdf5"):
        #     group_label = KMeans(n_clusters=cent_num, random_state=rs).fit_predict(total_data[:,[0,2]])
        #     h5f = h5py.File(foreground_path_ori + "/jkf_label.hdf5", "w")
        #     h5f["/group_label"] = group_label.astype(dtype=numpy.intc)
        #     h5f.close()
        # else:
        #     h5f = h5py.File(foreground_path_ori + "/jkf_label.hdf5", "r")
        #     group_label = h5f["/group_label"][()]
        #     h5f.close()


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


        print("Time: %.2f sec. %d galaxies"%(t2-t1, total_num))
    comm.Barrier()

    if rank > 0:
        h5f = h5py.File(stack_file_path, "r")
        total_data = h5f["/data"][()]
        group_label = h5f["/group_label"][()]
        h5f.close()
    cent_num = group_label.max() + 1

    # assign the source into the artificial exposures
    min_src_num = 100

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


elif cmd == "prepare_background":

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    cpus = comm.Get_size()

    backz_min = 0#float(argv[2])
    backz_max = 3#float(argv[3])

    if rank == 0:
        if os.path.exists(result_cata_path + "/background"):
            shutil.rmtree(result_cata_path + "/background")
        os.makedirs(result_cata_path + "/background")

        log_inform = "H0: %.4f. Omega_m0: %.4f"%(H0, omega_m0)
        print(log_inform)

    comm.Barrier()

    total_expos = []
    with open(fourier_avail_expo_path, "r") as f:
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

        # selection
        idx1 = data[:, nstar_idx] >= nstar_thresh
        idx2 = data[:, flux2_alt_idx] >= flux2_alt_thresh
        idx3 = numpy.abs(data[:, gf1_idx]) <= gf1_thresh
        idx4 = numpy.abs(data[:, gf2_idx]) <= gf2_thresh
        idx5_0 = data[:, Zp_idx] > 0
        idx5_1 = data[:, Zp_idx] >= backz_min
        idx5_2 = data[:, Zp_idx] <= backz_max
        idx5 = idx5_0 & idx5_1 & idx5_2
        # idx6 = numpy.sqrt(data[:, gf1_idx]**2 + data[:, gf2_idx]**2) <= gf_thresh
        # idx6 = numpy.abs(data[:,dist_idx]) <= dist_thresh

        idx = idx1 & idx2 & idx3 & idx4 & idx5 #& idx6

        src_num = idx.sum()

        if src_num > 0:
            src_data = data[idx]

            ra = src_data[:, ra_idx]
            dec = src_data[:, dec_idx]
            # print(ra.max(), ra.min())
            # one of the two survey areas spans ra=0
            idx = ra < 100
            ra[idx] += 360

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

            # G1, G2, N, U, V, RA, RA_radian, DEC, DEC_radian, COS(DEC), SIN(DEC), Z, Z_ERR, COMOVING DISTANCE
            dst_data = numpy.zeros((src_num, 14), dtype=numpy.float32)

            dst_data[:, 0] = src_data[:, mg1_idx]
            dst_data[:, 1] = src_data[:, mg2_idx]
            dst_data[:, 2] = src_data[:, mn_idx]
            # the signs of U & V are different with that in the paper
            dst_data[:, 3] = -src_data[:, mu_idx]
            dst_data[:, 4] = -src_data[:, mv_idx]
            
            dst_data[:, 5] = src_data[:, ra_idx]
            dst_data[:, 6] = src_data[:, ra_idx]*deg2rad
            dst_data[:, 7] = src_data[:, dec_idx]
            dst_data[:, 8] = src_data[:, dec_idx]*deg2rad
            dst_data[:, 9] = numpy.cos(src_data[:, dec_idx]*deg2rad)
            dst_data[:, 10] = numpy.sin(src_data[:, dec_idx]*deg2rad)
            dst_data[:, 11] = src_data[:, Zp_idx]
            # idx_z_select = dst_data[:,8] < 0.25
           # print(idx_z_select.sum(), src_num)
           #  dst_data[:,8][idx_z_select] = dst_data[:,8][idx_z_select] + 0.3
            # replace the photoZ by spectral Z
            # idxz = src_data[:, Zs_idx] > 0
            # dst_data[:, 8][idxz] = src_data[:, Zs_idx][idxz]

            # dst_data[:, 12] = src_data[:, 3]*2#(src_data[:, Z_B_MAX_idx] - src_data[:,Z_B_MIN_idx])/2
            dst_data[:, 13] = cosmos.comoving_distance(src_data[:, Zp_idx]).value*H0/100 # in unit of Mpc/h

            expo_dst_path = result_cata_path + "/background/%s" %expo_name
            expo_dst_path_pi2 = parent_path_pi2 + "/background_omg_%.3f/%s" %(omega_m0, expo_name)
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
    bin_scale = 100000000

    # for correlation calculation
    if rank == 0:
        # set up G bins for gamma_t calculation
        field_name = []
        with open(fourier_avail_expo_path, "r") as f:
            f_lines = f.readlines()
        for ff in f_lines:
            field_name.append(ff.split("\n")[0])

        # set up bins for G1(or G2) for PDF_SYM
        file_list = numpy.arange(0, len(field_name),dtype=numpy.intc)
        numpy.random.shuffle(file_list)

        for tag, i in enumerate(file_list[:1000]):

            # print(field_name[i])
            h5f_src = h5py.File(field_name[i], "r")
            temp = h5f_src["/data"][()][:, mg1_idx:mg2_idx + 1]
            if tag == 0:
                src_data = temp
            else:
                src_data = numpy.row_stack((src_data, temp))
            h5f_src.close()

        # G bins for tangential shear calculation
        mg_bin = hk_tool_box.set_bin(src_data[:, 0], mg_bin_num, bin_scale)



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
        fore_file_label = numpy.arange(0, foreground_num, dtype=numpy.intc)
        numpy.random.shuffle(fore_file_label)

        for fore_tag, fore_label in enumerate(fore_file_label):

            fore_c = fore_contents[fore_label]
            cc = fore_c.split("\n")[0].split("\t")
            fore_expo_path = cc[0]
            h5f = h5py.File(fore_expo_path, "r")
            fore_data = h5f["/data"][()]
            h5f.close()

            fore_c_num = min(5, fore_data.shape[0])

            for ic in range(fore_c_num):
                if Gt_num >= Gts_total_num:
                    break
                ic_z = fore_data[ic,6] # redshift
                ic_com_dist = fore_data[ic, 7] # com_dist
                ic_skypos = SkyCoord(ra=fore_data[ic,0]*units.deg, dec=fore_data[ic,2]*units.deg, frame="fk5")
                separation_angle = ic_skypos.separation(back_expo_skypos).deg
                idx = separation_angle < 1
                if idx.sum() > 0:
                    for ib in back_expo_labels[idx]:

                        if Gt_num >= Gts_total_num:
                            break

                        h5f = h5py.File(back_expo_path[ib], "r")
                        back_data = h5f["/data"][()]
                        h5f.close()

                        idxz = back_data[:,11] > ic_z + 0.1
                        sub_num = idxz.sum()
                        if sub_num > 0:
                            sub_data = back_data[idxz]
                            ib_skypos = SkyCoord(ra=sub_data[:,5]*units.deg, dec=sub_data[:,7]*units.deg,frame="fk5")
                            position_ang = 2*ic_skypos.position_angle(ib_skypos).radian
                            mgt = sub_data[:,0]*numpy.cos(position_ang) - sub_data[:,1]*numpy.sin(position_ang)
                            mnu = sub_data[:,2] + sub_data[:,3]*numpy.cos(2*position_ang) - sub_data[:,4]*numpy.sin(2*position_ang)

                            # print(ic_com_dist, sub_data[:, 10] - ic_com_dist)
                            # 1662895.2007121066 = c^2/4/pi/G  [M_sum/pc] /10^6 [h/pc]
                            # idx_zero = sub_data[:,10]-ic_com_dist == 0

                            # print(idx_zero.sum(),ic_z, sub_data[:,10][idx_zero],sub_data[:,8][idx_zero])
                            mgt = mgt*sub_data[:,13]/ic_com_dist/(sub_data[:,13]-ic_com_dist)*1662895.2007121066

                            if Gt_num + sub_num > Gts_total_num:
                                Gts[Gt_num: Gts_total_num] = mgt[Gts_total_num-Gt_num]
                                NU[Gt_num: Gts_total_num] = mnu[Gts_total_num-Gt_num]
                                Gt_num = Gts_total_num
                                # print("Finish: ",Gt_num)
                                # print(mgt[Gts_total_num-Gt_num][:100],mgt[Gts_total_num-Gt_num][-100:])
                                break
                            else:
                                Gts[Gt_num: Gt_num + sub_num] = mgt
                                NU[Gt_num: Gt_num + sub_num] = mnu
                                Gt_num += sub_num
                                # print("Finding: ",Gt_num)

            # print(fore_tag, Gt_num)
        mg_sigma_bin = hk_tool_box.set_bin(Gts, mg_bin_num, bin_scale)
        # print(tan_shear_guess)
        # print(delta_sigma_guess)
        print(mg_bin)
        print(mg_sigma_bin)
        # print(Gts)
        # print(Gts[:100],Gts[-100:])
        # print(numpy.abs(Gts).min())
        # h5f = h5py.File(result_cata_path + "/data4bin.hdf5", "w")
        # h5f["/Gt_sigma"] = Gts
        # h5f["/NU"] = NU
        # h5f["/Gt"] = src_data[:, 0]
        # h5f.close()

        # hist2d_mg_bin = numpy.zeros((hist2d_mg_num+1,))
        # hist2d_mnu_bin = numpy.zeros((hist2d_mg_num+1,))
        #
        # hist2d_mg_bin[:hist2d_mg_num2] = -hk_tool_box.set_bin_log(0.001, Gts.max()*500, hist2d_mg_num2)
        # hist2d_mg_bin[hist2d_mg_num2+1:] = hk_tool_box.set_bin_log(0.001, Gts.max()*500, hist2d_mg_num2)
        #
        # hist2d_mnu_bin[:hist2d_mg_num2] = -hk_tool_box.set_bin_log(0.001, NU.max()*500, hist2d_mg_num2)
        # hist2d_mnu_bin[hist2d_mg_num2+1:] = hk_tool_box.set_bin_log(0.001, NU.max()*500, hist2d_mg_num2)
        #
        # hist2d_mg_bin = numpy.sort(hist2d_mg_bin)
        # hist2d_mnu_bin = numpy.sort(hist2d_mnu_bin)
        #
        # hist2d_mg_bin[0] = hist2d_mg_bin[0]*bin_scale
        # hist2d_mg_bin[-1] = hist2d_mg_bin[-1]*bin_scale
        # hist2d_mnu_bin[0] = hist2d_mnu_bin[0]*bin_scale
        # hist2d_mnu_bin[-1] = hist2d_mnu_bin[-1]*bin_scale


        h5f = h5py.File(result_cata_path + "/pdf_inform.hdf5", "w")

        # h5f["/hist2d_mg_sigma_bin"] = hist2d_mg_bin.astype(dtype=numpy.float32)
        # h5f["/hist2d_mn_sigma_bin"] = hist2d_mnu_bin.astype(dtype=numpy.float32)

        h5f["/mg_sigma_bin"] = mg_sigma_bin.astype(dtype=numpy.float32)
        h5f["/mg_gt_bin"] = mg_bin.astype(dtype=numpy.float32)
        h5f["/gt_guess"] = tan_shear_guess
        h5f["/delta_sigma_guess"] = delta_sigma_guess
        h5f["/separation_bin"] = separation_bin[:19].astype(dtype=numpy.float32)
        h5f["/cosmological_params"] = numpy.array([H0, omega_m0], dtype=numpy.float32)

        h5f.close()

    comm.Barrier()
