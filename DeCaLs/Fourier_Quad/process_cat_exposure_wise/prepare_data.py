import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from mpi4py import MPI
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

if rank == 0:
    files = os.listdir("./g")
    target = []
    for fn in files:
        if ".hdf5" in fn:
            target.append(fn)
    file_list = tool_box.alloc(target, cpus)
    # file_list = [["%s"%j for i in range(j+1)] for j in range(cpus)]
    # print(file_list)
    print(len(target)," files")
else:
    file_list = None

sub_list = comm.scatter(file_list, root=0)
# print(sub_list)
print("%d gets %d files"%(rank, len(sub_list)))

for fn in sub_list:
    h5f = h5py.File("g/%s"%fn,"r")
    src = h5f["/data"][()]
    h5f.close()

    idx1 = src[:,2] < 24
    idx2 = src[:,6] < 0.05
    idx3 = src[:,14] >= 20
    idx4 = src[:,15] < 48
    idx5 = src[:,16] < 48
    idx6 = src[:,21] > 3
    idx7 = numpy.abs(src[:,24]) <= 0.0012
    idx8 = numpy.abs(src[:,25]) <= 0.0012

    idx = idx1 & idx2 & idx3 & idx4 & idx5 & idx6 & idx7 & idx8

    num = idx.sum()

    if num > 0:
        src_ = src[idx]
        dst = numpy.zeros((num, 8), dtype=numpy.float32)

        dst[:, 0] = src_[:,22]
        dst[:, 1] = src_[:,23]
        dst[:, 2] = src_[:,5]
        dst[:, 3:] = src_[:,26:]

        h5f = h5py.File("g_select/%s"%fn, "w")
        h5f["/data"] = dst
        h5f.close()


