import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
import h5py
from mpi4py import MPI
import matplotlib.pyplot as plt
import numpy
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


area_num = 4

result_source = "result_ext"

envs_path = "%s/work/envs/envs.dat"%my_home

gets_item = [["gg_lensing", "ggl_path_data", "0"]]
path_items = tool_box.config(envs_path, ["get"], gets_item)

data_path = path_items[0]

h5f_path_cut = data_path + "cata_%s_cut.hdf5" % result_source
h5f_path_grid = data_path + "cata_%s_grid.hdf5" % result_source

gets_item = [["fresh_para_idx", "nstar", "0"], ["fresh_para_idx", "flux_alt", "0"],
             ["fresh_para_idx", "ra", "0"], ["fresh_para_idx", "dec", "0"],
             ["fresh_para_idx", "gf1", "0"], ["fresh_para_idx", "gf2", "0"],
             ["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
             ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
             ["fresh_para_idx", "h2", "0"], ["fresh_para_idx", "total_area", "0"]]
gets = ["get" for i in range(len(gets_item))]
para_items = tool_box.config(envs_path, gets, gets_item)

nstar_lb = int(para_items[0])
flux_alt_lb = int(para_items[1])
total_area_lb = int(para_items[11])

ra_lb = int(para_items[2])
dec_lb = int(para_items[3])

field_g1_lb = int(para_items[4])
field_g2_lb = int(para_items[5])

mg1_lb = int(para_items[6])
mg2_lb = int(para_items[7])
mn_lb = int(para_items[8])
mu_lb = int(para_items[9])
mv_lb = int(para_items[10])
z_lb = -3


pre_h5f = h5py.File(h5f_path_cut,"r")
grid_h5f = h5py.File(h5f_path_grid,"r")


data_name = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V",
             "num_in_block", "block_start", "block_end", "block_boundy", "block_boundx"]


if rank < area_num:

    grid_data_sets = ["/background/w_%d/%s" % (rank + 1, data_name[i]) for i in range(len(data_name))]
    for i in range(8):
        pre_set_name = "/w_%d/%s"%(rank+1, data_name[i])
        grid_set_name = "/background/w_%d/%s"%(rank+1, data_name[i])
        attrs_name = "shape"
        pre_data = pre_h5f[pre_set_name].value
        grid_data = grid_h5f[grid_set_name].value
        grid_attrs = grid_h5f[grid_set_name].attrs[attrs_name]

        sp1 = pre_data.shape
        sp2 = grid_data.shape
        if sp1[0] != sp2[0]:
            print("The grid data doesn't match the original one")
        if grid_attrs[0] != sp1[0] or grid_attrs[0] != sp2[0]:
            print("The grid data shape doesn't match the original one")

comm.Barrier()

for area_id in range(1, area_num+1):

    pre_set_names = ["/w_%d/%s" % (area_id, data_name[i]) for i in range(8)]
    grid_set_names = ["/background/w_%d/%s" % (area_id, data_name[i])for i in range(len(data_name))]

    pre_datas = [pre_h5f[name].value for name in pre_set_names]
    grid_datas = [grid_h5f[grid_set_names[i]].value for i in range(8)]

    num_in_block = grid_h5f["/background/w_%d/%s"%(area_id, data_name[8])].value
    block_start = grid_h5f["/background/w_%d/%s"%(area_id, data_name[9])].value
    block_end = grid_h5f["/background/w_%d/%s" % (area_id, data_name[10])].value
    block_boundy = grid_h5f["/background/w_%d/%s" % (area_id, data_name[11])].value
    block_boundx = grid_h5f["/background/w_%d/%s" % (area_id, data_name[12])].value

    grid_shape = grid_h5f["/background/w_%d"%area_id].attrs["grid_shape"]
    grid_ny = grid_shape[0]
    grid_nx = grid_shape[1]
    grid_num = grid_ny * grid_nx
    tasks = [i for i in range(grid_num)]
    if rank == 0:
        print("Area: %d, grid: %d (%d, %d)"%(area_id,grid_num, grid_ny, grid_nx))
    my_grids = tool_box.allot(tasks, cpus)[rank]

    for ig in my_grids:
        row, col = divmod(ig, grid_nx)
        grid_s_id = block_start[ig, 0]
        grid_e_id = block_end[ig, 0]

        ra_min, ra_max = block_boundx[ig, 0], block_boundx[ig, 1]
        dec_min, dec_max = block_boundy[ig, 0], block_boundy[ig, 2]

        idx_1 = pre_datas[1] < ra_max
        idx_2 = pre_datas[1] >= ra_min
        idx_3 = pre_datas[2] < dec_max
        idx_4 = pre_datas[2] >= dec_min
        idx = idx_1 & idx_2 & idx_3 & idx_4
        for i in range(8):
            data_i = grid_datas[i][grid_s_id:grid_e_id]
            pre_data_i = pre_datas[i][idx]
            residual = numpy.sum(data_i-pre_data_i)
            if abs(residual) > 0.000001:
                print("Area: %d, %s grid[%d, %d] is wrong, residual: %f!!!"%(area_id, data_name[i], row, col, residual))




