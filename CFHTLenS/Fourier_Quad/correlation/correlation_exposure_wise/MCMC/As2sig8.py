import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import tool_box
import correlation_function_tool as cf_tool
import numpy
import time
import h5py
from mpi4py import MPI


t1 = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

data_path, H0 = argv[1], float(argv[2])

file_name = os.path.basename(data_path).split(".")[0]
parent_path = os.path.dirname(data_path)

As_col = 0
omega_cm0_col = 1
omega_bm0_col = 2
h = H0/100

# h5f = h5py.File(data_path,"r")
# src_data = h5f["/data"][()]
# h5f.close()
src_data = numpy.load(data_path)["arr_0"]

src_num, cols = src_data.shape

task_list = [1 for i in range(src_num)]
sub_task_list = tool_box.alloc(task_list, numprocs)
counts = [sum(sub_task_list[i]) for i in range(numprocs)]
sub_src_num = counts[rank]
row_st = [sum(counts[:i]) for i in range(numprocs)]
row_ed = [row_st[i] + counts[i] for i in range(numprocs)]

irow_st, irow_ed = row_st[rank], row_ed[rank]

# data_tran = src_data[irow_st:irow_ed]

cols = 3
data_tran = numpy.zeros((sub_src_num,cols))
data_tran[:,As_col] = src_data[irow_st:irow_ed,0]/10**9
data_tran[:,omega_cm0_col] = src_data[irow_st:irow_ed,1]*(1-src_data[irow_st:irow_ed,2])
data_tran[:,omega_bm0_col] = src_data[irow_st:irow_ed,1]*src_data[irow_st:irow_ed,2]

As = data_tran[:, As_col]
omega_cm0 = data_tran[:, omega_cm0_col]
omega_bm0 = data_tran[:, omega_bm0_col]
z = [0]

for i in range(numprocs):
    if rank == i:
        print("%d: %d-%d (%d) of %d"%(rank, irow_st, irow_ed, sub_src_num, src_num))
    comm.Barrier()
comm.Barrier()
sig8 = cf_tool.As2sigma8(As, omega_cm0, omega_bm0,z, H0)
data_tran[:, As_col] = sig8

comm.Barrier()

if rank > 0:
    comm.Send([data_tran, MPI.DOUBLE], dest=0, tag=rank)
else:

    final_data = numpy.zeros((src_num,cols))
    final_data[irow_st:irow_ed] = data_tran

    for ir in range(1, numprocs):
        irow_st, irow_ed = row_st[ir], row_ed[ir]
        recv_buf = numpy.empty((counts[ir], cols), dtype=numpy.float64)
        comm.Recv(recv_buf, source=ir, tag=ir)

        final_data[irow_st:irow_ed] = recv_buf

    numpy.savez(parent_path + "/" + file_name + "_s8.npz", final_data)


    diff_cm = src_data[:,1]*(1-src_data[:,2]) - final_data[:,omega_cm0_col]
    diff_bm = src_data[:,1]*src_data[:,2] - final_data[:,omega_bm0_col]
    print("Test: ",diff_cm.min(), diff_cm.max())
    print("Test: ",diff_bm.min(), diff_bm.max())
t2 = time.time()

comm.Barrier()
if rank == 0:
    print("H0: %.2f. Time: %.2f sec"%(H0, t2-t1))