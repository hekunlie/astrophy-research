import os
my_home = os.popen("echo $HK_MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import shutil
import h5py
import numpy
from mpi4py import MPI
import hk_tool_box
from astropy.cosmology import FlatLambdaCDM



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

# cosmological parameters
omega_m0 = float(argv[2])
H0 = 67.5
cosmos = FlatLambdaCDM(H0, omega_m0)
deg2arcmin = 60
deg2rad = numpy.pi/180

# with open("/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/haojie_cata/file_list", "r") as f:
# # with open("/home/hklee/work/catalog/Haojie_cata/lumin_z_bin_red_planck2018/file_list", "r") as f:
#     folder_name = f.readlines()

comm.Barrier()

parent_path = "/home/hklee/work/catalog/Haojie_cata"
# parent_path_pi2 = "/home/hklee/work/DECALS/DECALS_v210729/gg_lensing/cata"
# parent_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/haojie_cata"
parent_path_pi2 = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/DECALS_v210729/gg_lensing/cata"

folder_name = ["z_0.1_0.3_absz_-22.0_-23.0\n"]

min_src_num = 50#int(argv[2])

for fnm in folder_name:
    total_path = "%s/lumin_z_bin_red_planck2018/%s"%(parent_path, fnm.split("\n")[0])
    file_name = argv[1] + "_jkf.hdf5"

    final_path = parent_path_pi2 + "/%s_%s"%(argv[1],fnm.split("\n")[0])
    result_cata_path = parent_path + "/segment/%s_%s"%(argv[1],fnm.split("\n")[0])
    if rank == 0:
        if os.path.exists(result_cata_path):
            shutil.rmtree(result_cata_path)
        os.mkdir(result_cata_path)

    comm.Barrier()

    h5f = h5py.File("%s/%s"%(total_path, file_name), "r")
    src_data = h5f["/data"][()]
    pix_ra_dec = h5f["/pix_ra_dec"][()]
    pixel_count = h5f["/pix_count"][()]
    src_pix_label = h5f["/data_pix_label"][()]
    group_label = h5f["/jkf_label_200"][()]
    h5f.close()
    cent_num = group_label.max() + 1

    # src_num = total_data.shape[0]
    src_num = src_data.shape[0]
    total_data = numpy.zeros((src_num,  8), dtype=numpy.float32)
    total_data[:, 0] = src_data[:, 1]   # RA
    total_data[:, 1] = src_data[:, 1]*deg2rad   # RA_radian
    total_data[:, 2] = src_data[:, 2]   # Dec
    total_data[:, 3] = src_data[:, 2]*deg2rad   # Dec_radian
    total_data[:, 6] = src_data[:, 3]   # z
    h5f.close()

    total_data[:, 4] = numpy.cos(total_data[:, 3]) # Cos(DEC)
    total_data[:, 5] = numpy.sin(total_data[:, 3]) # SIN(DEC)
    total_data[:, 7] = cosmos.comoving_distance(total_data[:, 6]).value * H0 / 100  # Distance


    # assign the source into the artificial exposures
    expos_avail_sub = []
    expos_count = 0

    group_list = [i for i in range(cent_num)]
    print(group_list, cent_num)

    sub_group_list = hk_tool_box.alloc(group_list, cpus)[rank]

    # divide the group into many exposures
    for group_tag in sub_group_list:

        # select individual group
        idx_group = group_label == group_tag
        sub_data = total_data[idx_group]

        ground_src_num = idx_group.sum()

        if ground_src_num <= min_src_num:
            expos_name = "%d-0" %group_tag
            expos_path = result_cata_path + "/%s.hdf5" % expos_name
            expos_path_pi2 = final_path + "/%s.hdf5" % expos_name
            h5f_expos = h5py.File(expos_path, "w")
            h5f_expos["/data"] = sub_data
            h5f_expos.close()

            expos_avail_sub.append("%s\t%s\t%d\t%d\n"
                                   % (expos_path_pi2, expos_name, ground_src_num, group_tag))
            expos_count += 1
        else:
            m, n = divmod(ground_src_num, min_src_num)
            nums_distrib = hk_tool_box.alloc([1 for i in range(ground_src_num)], m)
            nums = [sum(nums_distrib[i]) for i in range(m)]
            nums_st = [sum(nums[:i]) for i in range(m)]
            for count in range(m):
                expos_name = "%d-%d" % (group_tag, count)
                expos_path = result_cata_path + "/%s.hdf5" % expos_name
                expos_path_pi2 = final_path + "/%s.hdf5" % expos_name
                h5f_expos = h5py.File(expos_path, "w")
                h5f_expos["/data"] = sub_data[nums_st[count]: nums_st[count]+nums[count]]
                h5f_expos.close()

                expos_avail_sub.append("%s\t%s\t%d\t%d\n"
                                       % (expos_path_pi2, expos_name, nums[count], group_tag))
                expos_count += 1

    comm.Barrier()
    expos_count_list = comm.gather(expos_count, root=0)
    expos_avail_list = comm.gather(expos_avail_sub, root=0)

    if rank == 0:

        buffer_expo = []
        for fns in expos_avail_list:
            buffer_expo.extend(fns)

        with open(result_cata_path + "/foreground_source_list.dat", "w") as f:
            f.writelines(buffer_expo)

        log_inform = "%d exposures\n"%len(buffer_expo)
        # with open("kmeans_log.dat", "w") as f:
        #     f.writelines(log_inform)
        print(result_cata_path, cent_num)
        print(log_inform)
        print(expos_count_list)
        print(sum(expos_count_list))

    comm.Barrier()
comm.Barrier()