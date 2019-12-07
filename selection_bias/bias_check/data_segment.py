import numpy
import h5py
from mpi4py import MPI
from sys import path, argv
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
import tool_box
import shutil
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

shear_num = 11
n,m = divmod(shear_num,numprocs)
tasks = [i for i in range(shear_num)]

my_task = tool_box.allot(tasks,numprocs)[rank]

print(rank, my_task)
data_type = argv[1]
sub_pts = 10
sub_num = 10000000
entry = [(sub_num*i, sub_num*(i+1)) for i in range(sub_pts)]


total_path = "/mnt/ddnfs/data_users/hkli/bias_check/data_from_pi/data_flux_pdf"
src_tag = 6
src_path = total_path + "/%d"%src_tag

src_shear_path = total_path + "/shear.hdf5"

for ig in my_task:

    data_path = total_path + "/data_%d_%s.hdf5"%(ig,data_type)
    h5f = h5py.File(data_path, "r")
    total_data_ig = h5f["/data"].value
    h5f.close()
    print("Reading total data of shear %d"%ig,total_data_ig.shape)

    for i in range(sub_pts):
        dst_path = total_path + "/sub_sample/%d/data_%d_%s.hdf5"%(i, ig, data_type)
        temp = total_data_ig[entry[i][0]:entry[i][1]]
        h5f = h5py.File(dst_path,"w")
        h5f["/data"] = temp
        h5f.close()

        if rank == 0:
            dst_shear_path = total_path + "/sub_sample/%d/shear.hdf5"%i
            if not os.path.exists(dst_shear_path):
                shutil.copyfile(src_shear_path, dst_shear_path)

        print("Write ",temp.shape," to %s"%dst_path)
