import numpy
import h5py
from mpi4py import MPI
from sys import path, argv
import shutil
import os

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

shear_num = 10
n,m = divmod(shear_num,numprocs)
tasks = [i for i in range(shear_num)]

my_task = [rank]

print(rank, my_task)
data_type = ["cross_term", "cross_term_est", "cross_term_sqrt",
             "noise_free","noise_residual","noisy_cpp"]

sub_pts = 10
sub_num = 20000000
entry = [(sub_num*i, sub_num*(i+1)) for i in range(sub_pts)]


total_path = argv[1]


h5f = h5py.File(total_path+"/data/shear.hdf5","r")
g1 = h5f["/g1"][()]
g2 = h5f["/g2"][()]
h5f.close()

if rank == 0:
    h5f = h5py.File(total_path + "/data_float/shear.hdf5", "w")
    h5f.create_dataset("/g1", data=g1, dtype=numpy.float32)
    h5f.create_dataset("/g2", data=g2, dtype=numpy.float32)
    h5f.close()

    for i in range(10):
        h5f = h5py.File(total_path + "/data_float/%d/shear.hdf5"%i, "w")
        h5f.create_dataset("/g1", data=g1, dtype=numpy.float32)
        h5f.create_dataset("/g2", data=g2, dtype=numpy.float32)
        h5f.close()


for ig in my_task:

    for dt in data_type:

        if "epsf" in total_path:
            data_path = total_path + "/data/data_%s_epsf_%d.hdf5"%(dt, ig)
            src_data_path = total_path + "/data_float/data_%s_epsf_%d.hdf5"%(dt, ig)
        else:
            data_path = total_path + "/data/data_%s_%d.hdf5"%(dt, ig)
            src_data_path = total_path + "/data_float/data_%s_%d.hdf5"%(dt, ig)

        h5f = h5py.File(data_path, "r")
        total_data_ig = h5f["/data"][()]
        h5f.close()
        print("Reading total data of shear %d"%ig, total_data_ig.shape)

        h5f = h5py.File(src_data_path, "w")
        h5f.create_dataset("/data", data=total_data_ig, dtype=numpy.float32)
        h5f.close()


        for i in range(sub_pts):

            if "epsf" in total_path:
                dst_path = total_path + "/data_float/%d/data_%s_epsf_%d.hdf5"%(i,dt,ig)
            else:
                dst_path = total_path + "/data_float/%d/data_%s_%d.hdf5"%(i,dt,ig)

            h5f = h5py.File(dst_path,"w")
            h5f.create_dataset("/data",data=total_data_ig[entry[i][0]:entry[i][1]], dtype=numpy.float32)
            h5f.close()

            print("Write to %s"%dst_path)
