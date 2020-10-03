import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
import matplotlib
# matplotlib.use("Agg")
# from plot_tool import Image_Plot
from mpi4py import MPI
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

total_path = argv[1]
# throw away
# sub_field_tag = int(argv[2])

src_path = total_path + "/fourier_cata"

source_list_nm = "nname_field_raw_avail.dat"

ab_fields = ["w1p4m0","w1p3p2","w1p3p3","w1p3m0","w1p3m4","w1p4p3","w1p2p1","w3p2m3","w3p1p2"]

fields = []


with open(src_path + "/cat_inform/"+source_list_nm, "r") as f:
    conts = f.readlines()

for nm in conts:
    field_nm = nm.split("\n")[0]
    if field_nm not in fields and field_nm not in ab_fields:
        fields.append(src_path + "/%s/result/%s_raw.hdf5"%(field_nm,field_nm))

if rank == 0:
    print("Totally %d fields" % len(fields))

# my_sub_area_list = []
# fnms = os.listdir(src_path + "/%s/result"%sub_fields[0])
# for nm in fnms:
#     if "p_all_raw" in nm:
#         my_sub_area_list.append(src_path + "/%s/result/%s"%(sub_fields[0],nm))

# if rank == 0:
#     print("Totally %d fields" % len(fields))
#     print(sub_field_tag,len(sub_fields)," fields, ", len(my_sub_area_list)," exposures",len(conts), expo_count, sub_fields)
# exit()
my_sub_area_list = tool_box.alloc(fields, cpus)[rank]

if len(my_sub_area_list) > 0:
    for tag, expo_path in enumerate(my_sub_area_list):
        h5f = h5py.File(expo_path,"r")
        temp = h5f["/data"][()]

        if tag == 0:
            data = temp
        else:
            data = numpy.row_stack((data, temp))
        h5f.close()

    sp = data.shape
else:
    sp = (0,0)
# print(rank,len(my_sub_area_list), len(fields), sp)

sp_total = comm.gather(sp, root=0)

comm.Barrier()
# if rank == 0:
#     print(rank,sp_total)
# exit()
if rank > 0 and sp[0] > 0:
    comm.Send([data,MPI.FLOAT], dest=0, tag=rank)
# if rank > 0:
#     pass
else:
    for ir in range(1, cpus):
        if sp_total[ir][0] > 0:
            recv_buf = numpy.empty(sp_total[ir],dtype=numpy.float32)
            comm.Recv(recv_buf,source=ir, tag=ir)
            data = numpy.row_stack((data, recv_buf))

    # for tag, expo_path in enumerate(my_sub_area_list):
    #
    #     h5f = h5py.File(expo_path, "r")
    #     temp = h5f["/data"][()]
    #
    #     if tag == 0:
    #         data = temp
    #     else:
    #         data = numpy.row_stack((data, temp))
    #     h5f.close()

    # h5f = h5py.File(total_path + "/selection_bias/anamoly_data/data_%d.hdf5" % sub_field_tag, "w")
    # h5f["/data"] = data
    # h5f.close()

    col_shift = 0

    nstar = data[:,col_shift+5]
    imax = data[:,col_shift+6]
    jmax = data[:,col_shift+7]
    gf1 = data[:, col_shift+15]
    gf2 = data[:, col_shift+16]

    idx1 = nstar >= 12
    idx2 = imax < 48
    idx3 = jmax < 48
    idx4 = numpy.abs(gf1) <= 0.005
    idx5 = numpy.abs(gf2) <= 0.005
    idx_ = idx1 & idx2 & idx3 & idx4 & idx5

    src_num = idx_.sum()

    data_sub = data[idx_]

    # sigma = data_sub[:, col_shift+3]
    # hlr = data_sub[:,col_shift+7]
    # hla = data_sub[:,col_shift+8]
    # snr = hlr/numpy.sqrt(hla)
    # snr_sort = numpy.sort(snr)
    # # snr_cut = [snr_sort[int(i*0.1*src_num)] for i in range(10)]
    # snr_cut = [0, snr.max()]
    #
    # flux2 = data_sub[:, col_shift+10]
    # flux2_alt = data_sub[:, col_shift+11]
    # flux2_alt_sort = numpy.sort(flux2_alt)
    # flux2_alt_cut = [flux2_alt_sort[int(i*0.1*src_num)] for i in range(10)]
    # flux2_alt_cut = [flux2_alt_sort[int(i*0.1*src_num)] for i in range(10)]

    gf1 = data_sub[:, col_shift+15]
    gf2 = data_sub[:, col_shift+16]

    data_fq = data_sub[:,col_shift+17:]


    bin_num1 = 50
    bin_num2 = 50
    gf1_bin = numpy.linspace(-0.005, 0.005, bin_num1+1,dtype=numpy.float32)
    gf2_bin = numpy.linspace(-0.005, 0.005, bin_num2+1,dtype=numpy.float32)
    gf1_pts = (gf1_bin[1:] + gf1_bin[:-1])/2
    gf2_pts = (gf2_bin[1:] + gf2_bin[:-1])/2

    # h5f = h5py.File(total_path + "/selection_bias/data/cutoff_%d.hdf5"%sub_field_tag, "w")
    h5f = h5py.File(total_path + "/selection_bias/cutoff.hdf5", "w")
    h5f["/gf1_bin"] = gf1_pts
    h5f["/gf2_bin"] = gf2_pts
    # h5f["/snr"] = snr_cut
    # h5f["/flux2_alt"] = flux2_alt_cut

    # print(gf1_pts)
    # print(gf2_pts)
    # print(h5f["/snr"][()])
    # print(h5f["/flux2_alt"][()])

    for i in range(bin_num1):
        idx_11 = gf1 >= gf1_bin[i]
        idx_12 = gf1 < gf1_bin[i+1]
        idx = idx_11 & idx_12
        num1 = idx.sum()

        h5f["/%d/mg_1"%i] = data_fq[:,0][idx]
        h5f["/%d/mn_1"%i] = data_fq[:,2][idx]
        h5f["/%d/mu_1"%i] = -data_fq[:,3][idx]

        # print(i, gf1_bin[i], gf1_bin[i + 1], num1)

    for i in range(bin_num2):
        idx_11 = gf2 >= gf2_bin[i]
        idx_12 = gf2 < gf2_bin[i+1]
        idx = idx_11 & idx_12
        num2 = idx.sum()

        h5f["/%d/mg_2"%i] = data_fq[:,1][idx]
        h5f["/%d/mn_2"%i] = data_fq[:,2][idx]
        h5f["/%d/mu_2"%i] = -data_fq[:,3][idx]

        # print(i, gf2_bin[i], gf2_bin[i + 1], num2)

    h5f.close()
comm.Barrier()