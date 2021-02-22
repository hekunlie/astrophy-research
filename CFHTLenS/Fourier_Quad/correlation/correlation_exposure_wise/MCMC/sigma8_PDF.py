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

As_col = 0
sig8_col = 1
omega_cm0_col = 2
omega_bm0_col = 3
h = 0.6775
H0 = h*100
z = [0]

num = int(argv[1])

data = numpy.zeros((num, 4))


rng = numpy.random.RandomState(rank*100+1)
As = rng.uniform(0.01,5,num)/10**9
omega_m0 = rng.uniform(0.05,0.7,num)
omega_bm0_ratio = rng.uniform(0.05,0.5,num)
omega_cm0 = omega_m0*omega_bm0_ratio*h*h
omega_bm0 = omega_m0*(1 - omega_bm0_ratio)*h*h


data[:,As_col] = As
data[:,omega_cm0_col] = omega_cm0
data[:,omega_bm0_col] = omega_bm0


for i in range(num):
    try:
        data[i,sig8_col] = cf_tool.get_CambResult(H0, omega_cm0[i], omega_bm0[i], As[i], 0.965, [0], kmax=3)[1]
    except:
        data[i,sig8_col] = - 10


comm.Barrier()

if rank > 0:
    comm.Send([data, MPI.DOUBLE], dest=0, tag=rank)
else:
    total_num = int(num*numprocs)
    final_data = numpy.zeros((num, 4))
    final_data[0:num] = data

    for ir in range(1, numprocs):
        irow_st, irow_ed = int(ir*num), int((ir+1)*num)
        recv_buf = numpy.empty((num, 4), dtype=numpy.float64)
        comm.Recv(recv_buf, source=ir, tag=ir)

        final_data[irow_st:irow_ed] = recv_buf
    idx = final_data[:,sig8_col] > -1
    cache = final_data[idx]
    numpy.savez("./sigma8_pdf.npz", cache, final_data)

t2 = time.time()

if rank == 0:
    print("%.2f sec"%(t2-t1))
comm.Barrier()