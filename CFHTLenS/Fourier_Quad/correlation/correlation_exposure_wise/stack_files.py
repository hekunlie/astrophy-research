import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import tool_box
import numpy
from mpi4py import MPI
import h5py


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

target_path = argv[1]
dst_path = argv[2]
name_need = argv[3]
name_exclude = argv[4]


if rank == 0:
    expos = []
    fns = os.listdir(target_path)
    print(len(fns), " files")
    for fn in fns:
        if name_need in fn and name_exclude not in fn:
            expos.append(target_path + "/" + fn)
else:
    expos = None

expos = comm.bcast(expos, root=0)

expos_num = len(expos)
if rank == 0:
    print(target_path)
    print(name_need)
    print(dst_path)
    print(expos_num, " exposures")

if expos_num < 1:
    exit()
my_sub_area_list = tool_box.alloc(expos, cpus)[rank]

# print(rank, i, len(my_sub_area_list))

if len(my_sub_area_list) > 0:
    for tag, expo_path in enumerate(my_sub_area_list):

        h5f = h5py.File(expo_path, "r")
        temp = h5f["/data"][()]
        h5f.close()
        # Nan check
        idx = numpy.isnan(temp)
        if idx.sum() > 0:
            num = temp.shape[0]
            label = numpy.arange(0,num)
            idx = numpy.isnan(temp[:,-2])
            print(label[idx][:5])
            print("Find Nan ", expo_path)
        if tag == 0:
            stack_data = temp
        else:
            stack_data = numpy.row_stack((stack_data, temp))

    sp = stack_data.shape
else:
    sp = (0, 0)

sp_total = comm.gather(sp, root=0)
comm.Barrier()

if rank > 0 and sp[0] > 0:
    comm.Send([stack_data, MPI.FLOAT], dest=0, tag=rank)
else:
    for ir in range(1, cpus):
        if sp_total[ir][0] > 0:
            recv_buf = numpy.empty(sp_total[ir], dtype=numpy.float32)
            comm.Recv(recv_buf, source=ir, tag=ir)
            stack_data = numpy.row_stack((stack_data, recv_buf))

    h5f = h5py.File(dst_path + "/stack_data.hdf5", "w")
    h5f["/data"] = stack_data
    h5f.close()
    h5f = h5py.File(dst_path + "/stack_data_ra_dec.hdf5", "w")
    h5f["/data"] = stack_data[:,:2]
    h5f.close()
comm.Barrier()
