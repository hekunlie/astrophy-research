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
from sys import argv


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]
cmds = ["collect","grid"]
if cmd not in cmds:
    if rank == 0:
        print("parameter must be one of ", cmds)
    exit()

area_num = 4

grid_scale = numpy.array([5, 10, 25, 60]) # arcmin
corre_scale = [10**(0.2*i) for i in range(11)]

cpu_block_size = int(cpus/area_num)

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "=" in path:
        env_location, env_path = path.split("=")[0:2]
        if "cfht_data_path" == env_location:
            data_path = env_path
        if "cfht_res_path" == env_location:
            res_path = env_path

# cpus_0 read the fields catalog and distribute them to related cpus
# cpu_0 & block_0, cpu_1 & block_1 ...
# cpu_i & block_i will be in charge of area_i

# data collection
if cmd == "collect":
    dicts, fields = tool_box.field_dict(data_path+"nname.dat")
    for i in range(area_num):
        if rank == 0:
            field_tar = []
            for field in fields:
                if "w%d"%(i+1) in field:
                    field_tar.append(field)
            field_pool = tool_box.allot(field_tar, cpus)
        else:
            field_pool = None

        field_pool = comm.scatter(field_pool, root=0)

        field_count = 0
        for field_name in field_pool:
            field_cat_path = data_path + "%s/result/%s_shear.dat"%(field_name, field_name)
            if os.path.exists(field_cat_path):
                try:
                    cat_arr = numpy.loadtxt(field_cat_path)
                    if field_count == 0:
                        cat_data = cat_arr
                    else:
                        cat_data = numpy.row_stack((cat_data, cat_arr))
                    field_count += 1
                except:
                    print(rank, "%s.dat doesn't exist"%field_name)
        data_sp = cat_data.shape
        data_sps = comm.gather(data_sp, root=0)

        if rank == 0:
            data_sps = numpy.array(data_sps)
            rows, cols = data_sps[:, 0], data_sps[:, 1]
            displ = []
            count = rows*cols
            for j in range(cpus):
                displ.append(count[0:j].sum())
            count = count.tolist()
            recv_buffer = numpy.empty((rows.sum(), cols[0]))
        else:
            count = None
            displ = None
            recv_buffer = None
        count = comm.bcast(count, root=0)
        displ = comm.bcast(displ, root=0)
        comm.Gatherv(cat_data, [recv_buffer, count, displ, MPI.DOUBLE], root=0)
        if rank == 0:
            final_data_path = res_path + "w%d.npz"%i
            numpy.savez(final_data_path, recv_buffer)

if cmd == "grid":
    if rank < area_num:
        # test
        cat_path = res_path + "w%d.npz"%rank + 1
        cat_data = numpy.load(cat_path)["arr_0"]
        ra, dec = cat_data[:,12]*60, cat_data[:,13]*60
        ra_min, ra_max, dec_min, dec_max = ra.min(), ra.max(), dec.min(), dec.max()
        grid_h5_path = res_path + "w%d_grid.hdf5"%rank + 1
        f = h5py.File(grid_h5_path, "w")
        for i, scale in enumerate(grid_scale):
            rows, cols = int((ra_max - ra_min)/scale+1), int((dec_max - dec_min)/scale+1)
            for row in range(rows):
                for col in range(cols):
                    group_name = "/%d/%d/%d"%(scale, row, col)
                    f.create_group(group_name)
                    idx1 = ra > ra_min + i*scale
                    f[group_name].attrs[]









