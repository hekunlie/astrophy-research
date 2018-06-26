import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import numpy
from mpi4py import MPI
import time
import tool_box
import warnings

warnings.filterwarnings("error")

# stack the CFHT shear catalog files into one big .npz file.

# if the 'filtered.dat' exists (includes the fields that will be used),
# it will stack the data of those fields in the .npz file.

# else, it will stack data from all the fields.
# it will produce the npz file of each field and each exposure of each field.

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()
with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        total_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_field_path" in path:
        field_path = path.split("=")[1]


data_cache = result_path + "data_cache.npz"

nname_path = total_path + "nname.dat"
field_dict, contents = tool_box.field_dict(nname_path)

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

fcount = 0
for i in range(len(path_list)):
    data_path = path_list[i]
    if os.path.exists(data_path):
        if filter_exist:
            temp = numpy.load(data_path)["arr_0"]
        else:
            area = os.path.basename(data_path)
            temp = numpy.loadtxt(data_path)
            field_name = area.split("_")[0]
            w_path = result_path + "field/%s/"%field_name
            if not os.path.exists(w_path):
                os.mkdir(w_path)
            area_npz = w_path + "%s.npz"%area
            npz_name = field_name + "\n"
            npz_cat.append(npz_name)
            numpy.savez(area_npz, temp)
            # for expo in list(field_dict[field_name].keys()):
            #     chip_count = 0
            #     for chip in field_dict[field_name][expo]:
            #         chip_shear_path = total_path + "%s/result/%s_shear.dat"%(field_name, chip)
            #         try:
            #             chip_data_temp = numpy.loadtxt(chip_shear_path,skiprows=1)
            #             if chip_count == 0:
            #                 expo_data = chip_data_temp.copy()
            #             else:
            #                 expo_data = numpy.row_stack((expo_data, chip_data_temp))
            #             chip_count += 1
            #         except:
            #             print("Empty %s/%s.shear.dat"%(field_name,chip))
            #     expo_data_path = field_path + "%s/%s.npz"%(field_name,expo)
            #     numpy.savez(expo_data_path, expo_data)
        if fcount == 0:
            data = temp.copy()
        else:
            data = numpy.row_stack((data, temp))
        fcount += 1
    else:
        print(rank, os.path.basename(data_path))

recv_fcount = comm.gather(fcount,root=0)
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
    print("Totally, %d galaxies are detected in %d (%d) fields"%(len(data), sum(recv_fcount), field_count))

te = time.clock()
if rank == 0:
    # if not filter_exist:
    #     npzs = []
    #     for sub_ in recv_npz_cat:
    #         npzs.extend(sub_)
    #     with open(filter_path, "w") as f:
    #         f.writelines(npzs)
    print(te-ts)