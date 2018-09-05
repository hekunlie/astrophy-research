import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
import h5py
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

area_num = 4
grid_scale = [1.6] # arcmin

cpu_block_size = int(cpus/area_num)

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "=" in path:
        env_location, env_path = path.split("=")[0:2]
        if "cfht_data_path" == env_location:
            data_path = env_path

# cpus_0 read the fields catalog and distribute them to related cpus
# cpu_0 & block_0, cpu_1 & block_1 ...
# cpu_i & block_i will be in charge of area_i
dicts, fields = tool_box.field_dict(data_path+"nname.dat")
for i in range(area_num):
    if rank == 0:
        field_tar = []
        for field in fields:
            if "w%d"%(i+1) in field:
                field_tar.append(field)
        field_pool = tool_box.allot(field_tar, cpus)

        print(" w%d field counts:"%(i+1), len(field_tar))
    else:
        field_pool = None

    field_pool = comm.scatter(field_pool, root=0)
    print(rank, field_pool, len(field_pool))

    field_count = 0
    for field_name in field_pool:
        field_cat_path = data_path + "%s/%s.dat"%(field_name, field_name)
        if os.path.exists(field_cat_path):
            try:
                cat_arr = numpy.loadtxt(field_cat_path)
                if field_count == 0:
                    cat_data = cat_arr
                else:
                    cat_data = numpy.row_stack((cat_data, ))
                field_count += 1
            except:
                print(rank, "%s.dat doesn't exist"%field_name)
    data_sp = cat_data.shape
    data_sps = comm.gather(data_sp, root=0)

    if rank == 0:
        rows = numpy.array(data_sps)[:,0].sum()
        recv_buffer = numpy.empty(())
#









