import matplotlib
matplotlib.use('Agg')
import os
from sys import path
path.append('/home/hkli/work/fourier_quad')
import numpy
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()
with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data" in path:
        total_path = path.split("=")[1]
    elif "cfht_res" in path:
        result_path = path.split("=")[1]


data_cache = result_path + "data_cache.npz"

nname_path = total_path + "nname.dat"
with open(nname_path,'r') as f:
    contents = f.readlines()

area_paths = []
for area in contents:
    if "w" in area:
        area = area.split("\n")[0]
        area_path = total_path + area + "/result/%s_shear.dat"%area
        area_paths.append(area_path)

m,n = divmod(len(area_paths), cpus)

if rank < cpus - 1:

    for i in range(m):
        if i == 0:
            data = numpy.loadtxt(area_paths[rank*m+i])
        else:
            data = numpy.row_stack((data, numpy.loadtxt(area_paths[rank*m+i])))
else:
    for i in range(m+n):

        if i == 0:
            data = numpy.loadtxt(area_paths[rank*m+i])
        else:
            data = numpy.row_stack((data, numpy.loadtxt(area_paths[rank*m+i])))
sp = data.shape
recv_sp = comm.gather(sp, root=0)

if rank > 0:
    comm.Send(data, dest=0, tag=rank)
else:
    for procs in range(1, cpus):
        recvs = numpy.empty(recv_sp[procs], dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        data = numpy.row_stack((data, recvs))

    numpy.savez(data_cache, data)
    print("Totally, %d galaxies are detected"%len(data))
te = time.clock()
if rank == 0 :
    print(te-ts)