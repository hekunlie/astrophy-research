import matplotlib
matplotlib.use('Agg')
import os
from sys import path
path.append('/home/hkli/work/fourier_quad')
import numpy
from mpi4py import MPI
import time
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()
with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        total_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]


data_cache = result_path + "data_cache.npz"

filter_path = result_path + "field/filtered.dat"
filter_exist = os.path.exists(filter_path)
if filter_exist:
    if rank==0:
        print("Stack filtered data")
    with open(filter_path, "r") as f:
        contents = f.readlines()
else:
    if rank==0:
        print("Stack unfiltered data")
    nname_path = total_path + "nname.dat"
    contents = tool_box.field_dict(nname_path)[1]

area_paths = []
field_count = 0
for area in contents:
    area = area.split("\n")[0]
    if filter_exist:
        area_path = result_path + "field/%s/%s_shear.dat.npz" %(area,area)
        field_count += 1
    else:
        area_path = total_path + area + "/result/%s_shear.dat"%area
        field_count += 1
    area_paths.append(area_path)

path_list = tool_box.allot(area_paths, cpus)[rank]
npz_cat = []

for i in range(len(path_list)):
    data_path = path_list[i]
    if filter_exist:
        temp = numpy.load(data_path)["arr_0"]
    else:
        area = os.path.basename(data_path)
        temp = numpy.loadtxt(data_path)
        w_path = result_path + "field/%s/"%area.split("_")[0]
        if not os.path.exists(w_path):
            os.mkdir(w_path)
        area_npz = w_path + "%s.npz"%area
        npz_name = area.split("_")[0] + "\n"
        npz_cat.append(npz_name)
        numpy.savez(area_npz, temp)
    if i == 0:
        data = temp
    else:
        data = numpy.row_stack((data, temp))

sp = data.shape
recv_sp = comm.gather(sp, root=0)
recv_npz_cat = comm.gather(npz_cat, root=0)

if rank > 0:
    comm.Send(data, dest=0, tag=rank)
else:
    for procs in range(1, cpus):
        recvs = numpy.empty(recv_sp[procs], dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        data = numpy.row_stack((data, recvs))
    numpy.savez(data_cache, data)
    print("Totally, %d galaxies are detected in %d fields"%(len(data), field_count))

te = time.clock()
if rank == 0:
    if not filter_exist:
        npzs = []
        for sub_ in recv_npz_cat:
            npzs.extend(sub_)
        with open(filter_path, "w") as f:
            f.writelines(npzs)
    print(te-ts)