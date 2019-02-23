import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
import h5py
from mpi4py import MPI
from sys import argv
import matplotlib.pyplot as plt
import numpy
import time

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# cpus = comm.Get_size()

area_num = 4

result_source = "result_int"
sect = "cfht"
envs_path = "%s/work/envs/envs.dat"%my_home

gets_item = [["cfht", "cfht_path_data", "0"]]
path_items = tool_box.config(envs_path, ["get"], gets_item)

data_path = path_items[0]
cf_cata_data_path = data_path + "cf_cata_%s_multi_.hdf5" % result_source

h5f = h5py.File(cf_cata_data_path, "r")
logger = tool_box.get_logger("./check_log.dat")
for area_id in range(area_num):

    grid_shape = h5f["/grid/w_%d"%area_id].attrs["grid_shape"]
    grid_ny, grid_nx = grid_shape
    grid_num = grid_ny*grid_nx
    print("\n")
    print(grid_shape)
    num_in_block = 0
    block_start = 0
    block_end = 0
    set_nms = ["boundy", "boundx", "num_in_block","block_start","block_end","G1_bin", "G1_bin"]
    for set_name in set_nms:
        data = h5f["/grid/w_%d/%s"%(area_id,set_name)].value
        data_sp = data.shape
        sp = h5f["/grid/w_%d/%s"%(area_id,set_name)].attrs["shape"]
        if set_name == "num_in_block":
            num_in_block = data
        if set_name == "block_start":
            block_start = data
        if set_name == "block_end":
            block_end = data
        if len(sp) == 1:
            log_info = "%s: data_shape: [%d, ], shape from attrs: [%d, ]"%(set_name, data_sp[0], sp[0])
        else:
            log_info = "%s: data_shape: [%d, %d], shape from attrs: [%d, %d]" \
                       % (set_name, data_sp[0], data_sp[1], sp[0], sp[1])
        logger.info(log_info)
        print(set_name, data_sp, sp)

    set_nms = ["data/G1","data/G2","data/N","data/U","data/V","data/RA","data/DEC"]
    for set_name in set_nms:
        data = h5f["/grid/w_%d/%s"%(area_id,set_name)].value
        data_sp = data.shape
        sp = h5f["/grid/w_%d/%s"%(area_id,set_name)].attrs["shape"]
        if len(sp) == 1:
            log_info = "%s: data_shape: [%d, ], shape from attrs: [%d, ]"%(set_name, data_sp[0], sp[0])
        else:
            log_info = "%s: data_shape: [%d, %d], shape from attrs: [%d, %d]" \
                       % (set_name, data_sp[0], data_sp[1], sp[0], sp[1])

        logger.info(log_info)
        log_info = "%s, num_check: num_in_block - data_shape: %d; data_shape - attrs_shape: %d"\
                   %(set_name, numpy.sum(num_in_block)-data_sp[0], data_sp[0]-sp[0])
        print(set_name, numpy.sum(num_in_block)-data_sp[0], data_sp[0]-sp[0])

    check_num = numpy.zeros((grid_num, 2), dtype=numpy.intc)
    print(list(h5f["/grid/w_%d/"%area_id].keys()))
    for row in range(grid_ny):
        for col in range(grid_nx):
            grid_id = row*grid_nx + col
            set_name = "/grid/w_%d/row_%d/col_%d"%(area_id,row, col)
            print(set_name)
            if num_in_block[grid_id] != 0:
                data = h5f[set_name].value
                data_sp = data.shape
            else:
                data_sp = (0,0)
            sp = h5f[set_name].attrs["shape"]
            check_num[grid_id] = data_sp[0], sp[0]
    print("data number check: %d, %d"%(numpy.sum(check_num[:,0] - num_in_block), numpy.sum(check_num[:,1] - num_in_block)))