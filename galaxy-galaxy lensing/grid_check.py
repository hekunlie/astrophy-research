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

logger = tool_box.get_logger("/home/hklee/work/CFHT/gg_lensing/log/check_%d.dat"%rank)

pre_h5f = h5py.File(h5f_path_cut, "r")
grid_h5f = h5py.File(h5f_path_grid, "r")


data_name = ["Z", "DISTANCE", "RA", "DEC", "COS_DEC", "G1", "G2", "N", "U", "V",
             "num_in_block", "block_start", "block_end", "block_boundy", "block_boundx"]


if rank < area_num:

    grid_data_sets = ["/background/w_%d/%s" % (rank + 1, data_name[i]) for i in range(len(data_name))]
    for i in range(10):
        pre_set_name = "/w_%d/%s"%(rank+1, data_name[i])
        grid_set_name = "/background/w_%d/%s"%(rank+1, data_name[i])
        attrs_name = "shape"
        pre_data = pre_h5f[pre_set_name].value
        grid_data = grid_h5f[grid_set_name].value
        grid_attrs = grid_h5f[grid_set_name].attrs[attrs_name]

        sp1 = pre_data.shape
        sp2 = grid_data.shape
        if sp1[0] != sp2[0]:
            logger.info("The grid data doesn't match the original one")
            print("The grid data doesn't match the original one")
        if grid_attrs[0] != sp1[0] or grid_attrs[0] != sp2[0]:
            logger.info("The grid data shape doesn't match the original one")
            print("The grid data shape doesn't match the original one")

comm.Barrier()

for area_id in range(1, area_num+1):

    pre_set_names = ["/w_%d/%s" % (area_id, data_name[i]) for i in range(10)]
    grid_set_names = ["/background/w_%d/%s" % (area_id, data_name[i])for i in range(len(data_name))]

    pre_RA = pre_h5f["/w_%d/RA"%area_id].value
    pre_DEC = pre_h5f["/w_%d/DEC" % area_id].value

    num_in_block = grid_h5f["/background/w_%d/num_in_block"%area_id].value
    block_start = grid_h5f["/background/w_%d/block_start"%area_id].value
    block_end = grid_h5f["/background/w_%d/block_end"%area_id].value
    block_boundy = grid_h5f["/background/w_%d/block_boundy"%area_id].value
    block_boundx = grid_h5f["/background/w_%d/block_boundx"%area_id].value

    dec_bin = grid_h5f["/background/w_%d/DEC_bin"%area_id].value
    ra_bin = grid_h5f["/background/w_%d/RA_bin"%area_id].value

    grid_shape = grid_h5f["/background/w_%d"%area_id].attrs["grid_shape"]
    grid_ny = grid_shape[0]
    grid_nx = grid_shape[1]
    grid_num = grid_ny * grid_nx
    tasks = [i for i in range(grid_num)]
    if rank == 0:
        print("Area: %d, grid: %d (%d, %d)"%(area_id,grid_num, grid_ny, grid_nx))
        logger.info("Area: %d, grid: %d (%d, %d)"%(area_id,grid_num, grid_ny, grid_nx))
    my_grids = tool_box.allot(tasks, cpus)[rank]

    # check col by col
    for i in range(10):
        # the original data and the grid data
        pre_datas = pre_h5f[pre_set_names[i]].value
        grid_datas = grid_h5f[grid_set_names[i]].value

        logger.info("Area: %d, %s" % (area_id, data_name[i]))

        for ig in my_grids:
            row, col = divmod(ig, grid_nx)
            block_s_id = block_start[ig, 0]
            block_e_id = block_end[ig, 0]

            grod_data_i = grid_datas[block_s_id:block_e_id]

            ra_min_b, ra_max_b = block_boundx[ig, 0], block_boundx[ig, 1]
            dec_min_b, dec_max_b = block_boundy[ig, 0], block_boundy[ig, 2]

            ra_min, ra_max = ra_bin[col,0], ra_bin[col+1,0]
            dec_min, dec_max = dec_bin[row,0], dec_bin[row + 1,0]

            check_1 = sum([abs(ra_min_b - ra_min), abs(ra_max_b - ra_max),
                           abs(dec_min_b - dec_min), abs(dec_max_b - dec_max)])
            if check_1 > 0.00001:
                print("Area: %d, The block boundary and bins of RA and DEC are not consistent"%area_id)
                logger.info("Area: %d, The block boundary and bins of RA and DEC are not consistent"%area_id)
            idx_1 = pre_RA < ra_max
            idx_2 = pre_RA >= ra_min
            idx_3 = pre_DEC < dec_max
            idx_4 = pre_DEC >= dec_min
            idx = idx_1 & idx_2 & idx_3 & idx_4

            delta_i = numpy.sum(pre_datas[idx] - grod_data_i)

            if abs(delta_i) > 0.00001:
                print("Area: %d, %s in grid[%d, %d] is wrong, residual: %f!!!"%(area_id, data_name[i], row, col, delta_i))
                logger.info("Area: %d, %s in grid[%d, %d] is wrong, residual: %f!!!"%(area_id, data_name[i], row, col, delta_i))
        comm.Barrier()

    comm.Barrier()

pre_h5f.close()
grid_h5f.close()


