import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box
import warnings
import shutil


#
# total_path = argv[1]
#
# ori_path = "/lustre/home/acct-phyzj/phyzj/CFHT/i"
#
# files_nm = os.listdir(ori_path)
# field_nm = []
#
#
# for fnm in files_nm:
#     if "w" in fnm:
#
#         expo_files_nm = os.listdir(ori_path + "/%s/result"%fnm)
#
#         temp_str = fnm
#         expo_tag = 0
#         for enm in expo_files_nm:
#             if "_all.cat" in enm:
#                 temp_str += "\t%s"%enm.split(".")[0]
#                 expo_tag += 1
#         temp_str += "\n"
#
#         field_nm.append(temp_str)
#
# with open("field_expo.dat", "w") as f:
#     f.writelines(field_nm)

total_path = argv[1]
dst_path = total_path + "/cat_hdf5/stack_expo"

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

with open("/home/hklee/work/CFHT/field_expo.dat","r") as f:
    field_expo = f.readlines()

if rank == 0:
    print(len(field_expo))

field_expo_sub = tool_box.alloc(field_expo,numprocs)[rank]

if rank == 0:
    if os.path.exists(dst_path):
        shutil.rmtree(dst_path,ignore_errors=True)
    os.mkdir(dst_path)

field_avail_sub = []

for fnm in field_expo_sub:
    temp_str = fnm.split()
    field_name = temp_str[0]

    max_len = 0
    data_list = []

    for tag, nm in enumerate(temp_str[1:]):
        expo_path = "%s/cat_hdf5/%s.hdf5" % (total_path, nm)

        if os.path.exists(expo_path):
            h5f = h5py.File(expo_path, "r")
            data = h5f["/data"][()]
            h5f.close()

            data_list.append(data)

            row, col = data.shape

            # the galaxy label in the field
            max_len = max(max_len, int(data[:, 17].max()))

            # print(nm, max_len)
        else:
            print(expo_path, os.path.exists(expo_path))

    if max_len > 0:

        data_stack = numpy.zeros((max_len, col), dtype=numpy.float32)
        count = numpy.zeros((max_len,), dtype=numpy.intc)
        galaxy_label = numpy.arange(0, max_len)
        # print(max_len)

        for i in range(len(data_list)):
            row = data_list[i].shape[0]
            ig = data_list[i][:, 17].astype(dtype=numpy.intc) - 1
            for j in range(row):
                i_gal = ig[j]
                count[i_gal] += 1
                data_stack[i_gal] += data_list[i][j]

        idx = count > 0
        # print(count.min(), count.max())

        galaxy_label = galaxy_label[idx]
        data_stack = data_stack[idx]

        final_data = numpy.zeros((idx.sum(), col))
        for i in range(col):
            final_data[:, i] = data_stack[:, i] / count[idx]

        dst_file = dst_path + "/%s.hdf5"%field_name
        h5f = h5py.File(dst_file,"w")
        h5f["/data"] = final_data
        h5f.close()

        field_avail_sub.append(dst_file + "\n")
    else:
        print("Empty: %s"%field_name)
comm.Barrier()

total_field_avail = comm.gather(field_avail_sub, root=0)

if rank == 0:

    field_avail = []
    for fa in total_field_avail:
        field_avail.extend(fa)

    with open(total_path + "/cat_inform/field_avail.dat","w") as f:
        f.writelines(field_avail)

    print(len(field_avail)," fields")