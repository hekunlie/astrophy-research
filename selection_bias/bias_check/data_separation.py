import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append('%s/work/mylib/' % my_home)
import h5py
from mpi4py import MPI
import tool_box
import numpy

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

data_path = argv[1]
if rank == 0:
    all_files = os.listdir(data_path + "/data_ori")

    files = []
    for al in all_files:
        if "data" in al:
            files.append(al)
else:
    files = None

comm.Barrier()

files = comm.bcast(files, root=0)
my_files = tool_box.alloc(files, numprocs)[rank]

for mf in my_files:
    h5f = h5py.File(data_path + "/data_ori/%s"%mf,"r")
    data = h5f["/data"][()]
    h5f.close()

    h5f = h5py.File(data_path + "/data/%s"%mf, "w")
    h5f["/mg1"] = data[:,0]
    h5f["/mg2"] = data[:,1]
    h5f["/mn"] = data[:,2]
    h5f["/mu"] = data[:,3]
    h5f["/mv"] = data[:,4]
    h5f.close()

comm.Barrier()