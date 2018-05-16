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

npz_path = result_path + "field/npz_files.dat"
npz_exist = os.path.exists(npz_path)

filter_path = result_path + "field/filtered.dat"
filter_exist = os.path.exists(filter_path)
if filter_exist:
    if rank==0:
        print("Stack filtered data")
    with open(filter_path, "r") as f:
        contents = f.readlines()
else:
    nname_path = total_path + "nname.dat"
    with open(nname_path, 'r') as f:
        contents = f.readlines()

area_paths = []
field_count = 0
for area in contents:
    if "w" in area:
        area = area.split("\n")[0]
        if filter_exist or npz_exist:
            area_path = result_path + "field/%s_shear.dat.npz" %area
            field_count += 1
        else:
            area_path = total_path + area + "/result/%s_shear.dat"%area
            field_count += 1
        area_paths.append(area_path)

m, n = divmod(len(area_paths), cpus)
npz_cat = []
if rank < cpus - 1:
    for i in range(m):
        data_path = area_paths[rank * m + i]
        if filter_exist or npz_exist:
            temp = numpy.load(data_path)["arr_0"]
        else:
            area = os.path.basename(data_path)
            temp = numpy.loadtxt(data_path)
            area_npz = result_path + "field/%s.npz"%area
            npz_cat.append(area+".npz")
            numpy.savez(area_npz, temp)
        if i == 0:
            data = temp
        else:
            data = numpy.row_stack((data, temp))
else:
    for i in range(m+n):
        data_path = area_paths[rank*m+i]
        if filter_exist or npz_exist:
            temp = numpy.load(data_path)["arr_0"]
        else:
            area = os.path.basename(data_path)
            temp = numpy.loadtxt(data_path)
            area_npz = result_path + "field/%s.npz"%area
            npz_cat.append(area + ".npz")
            numpy.savez(area_npz, temp)
        if i == 0:
            data = temp
        else:
            data = numpy.row_stack((data, temp))

sp = data.shape
recv_sp = comm.gather(sp, root=0)
if rank > 0:
    comm.Send(data, dest=0, tag=rank)
else:
    if not npz_exist:
        with open(npz_path,"w") as npz_f:
            npz_f.writelines(["npz"])

    for procs in range(1, cpus):
        recvs = numpy.empty(recv_sp[procs], dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        data = numpy.row_stack((data, recvs))

    numpy.savez(data_cache, data)

    print("Totally %d fields" % field_count)
    print("Totally, %d galaxies are detected"%len(data))
te = time.clock()
if rank == 0:
    print(te-ts)