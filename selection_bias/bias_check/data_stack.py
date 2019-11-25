import numpy
import h5py
from mpi4py import MPI
from sys import path
path.append("/home/hklee/work/mylib")
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

shear_num = 20
n,m = divmod(shear_num,numprocs)
tasks = [i for i in range(shear_num)]

my_task = tool_box.allot(tasks,numprocs)[rank]

print(rank, my_task)

sub_num = 10000000
total_data = numpy.zeros((sub_num*5,8))
entry = [(0,sub_num),(sub_num,3*sub_num),(3*sub_num, 5*sub_num)]

for ig in my_task:
    for i in range(3):
        data_path = "/mnt/perc/hklee/bias_check/result/data_%d/data_%d.hdf5"%(i, ig)
        h5f = h5py.File(data_path,"r")
        temp = h5f["/data"].value
        sp = temp.shape
        m,n = entry[i]
        total_data[m:n] = temp
        h5f.close()
        print(rank, "Reading shear %d data %d, (%d, %d) and sign to total data[%d, %d], %f, %f"
              % (ig, i, sp[0], sp[1],m,n,temp[0,2]-total_data[m,2], temp[sp[0]-1,2]-total_data[n-1,2]))

    sp = total_data.shape
    print("Write total data shear %d data %d, (%d, %d)"%(ig, ig, sp[0],sp[1]))
    total_data_path = "/mnt/perc/hklee/bias_check/result/datas/0/data_%d.hdf5" % ig
    h5f = h5py.File(total_data_path, "w")
    h5f["/data"] = total_data
    h5f.close()

    for i in range(1, 6):
        sub_data_path = "/mnt/perc/hklee/bias_check/result/datas/%d/data_%d.hdf5" % (i, ig)
        i = i - 1
        sub_data = total_data[i*sub_num:(i+1)*sub_num]
        h5f = h5py.File(sub_data_path,"w")
        h5f["/data"] = sub_data
        h5f.close()
