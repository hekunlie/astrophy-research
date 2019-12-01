import h5py
from mpi4py import MPI
from sys import path,argv
path.append("/home/hklee/work/mylib")
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

shear_num = 7
n, m = divmod(shear_num, numprocs)
tasks = [i for i in range(shear_num)]

my_task = tool_box.allot(tasks,numprocs)[rank]

print(rank, my_task)

total_path = argv[1]
# # columns: [col_st, col_ed]
# col_st, col_ed = int(argv[2]), int(argv[3])
# # name of the sub-data file
# sep_data_nm = argv[4]

dst_nms = ["gauss_noise_1.hdf5", "gauss_noise_2.hdf5", "gauss_noise_residual.hdf5",
           "moffat_noise_1.hdf5", "moffat_noise_2.hdf5", "moffat_noise_residual.hdf5"]

for ig in my_task:
    src_path = total_path + "/data_%d.hdf5"%ig
    h5f = h5py.File(src_path, "r")
    total_data_ig = h5f["/data"].value
    h5f.close()
    print("Reading total data of shear %d"%ig,total_data_ig.shape)

    for i in range(len(dst_nms)):
        dst_path = total_path + "/data_%d_%s"%(ig,dst_nms[i])
        h5f = h5py.File(dst_path,"w")
        h5f["/data"] = total_data_ig[:,i*4:(i+1)*4]
        h5f.close()
        print("Write data %s of shear %d" % (dst_nms[i],ig), total_data_ig.shape)