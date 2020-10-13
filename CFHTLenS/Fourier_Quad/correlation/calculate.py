import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import time
from mpi4py import MPI
from Fourier_Quad import *
import tool_box
import numpy
import h5py
import warnings
from sys import argv
from astropy.io import fits
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

envs_path = "%s/work/envs/envs.dat"%my_home
data_path, result_path = tool_box.config(envs_path, ['get','get'], [["correlation", "cfht_data_path", "1"],
                                                                    ["correlation", "cfht_result_path", "1"]])


the data and directories have not been prepared !!!

areas = ['w1', 'w2', 'w3', 'w4']
scale_num = 15
# aresec, the side length of each square block
block_scale = 5
radius_scales = [0]
for i in range(scale_num):
    radius_scales.append(10**(0.15 * i))

# read the number of the max length of the blocks.
# the data of each block will be stored in a line
# of the shared array in the memory.

data_h5_path = data_path + ".../"!!!
h5_file = h5py.File(data_h5_path,"r")
# max_block_length is maximum of galaxy number in one block
# of all sky area.
max_block_length = h5_file["grid"].attrs["max_block_length"][0,0]
# the maximum of block number of all sky area
max_block_num = h5_file["grid"].attrs["max_grid_length"][0,0]
h5_file.close()

# alloc the memory for data
size_tb = MPI.DOUBLE.Get_size()
size_ti = MPI.INT.Get_size()
if rank == 0:
    nb1 = int(7*max_block_num*max_block_length*size_tb)
    nb2 = int(max_block_num*max_block_length*size_tb)
else:
    nb1 = 0
    nb2 = 0
win_1 = MPI.Win.Allocate_shared(nb1, size_tb, comm=comm)
data_buf, size_tb = win_1.Shared_query(0)
win_2 = MPI.Win.Allocate_shared(nb2, size_ti, comm=comm)
label_buf, size_ti = win_2.Shared_query(0)
# each block occupies 7 contiguous rows
# G1, G2, N, U, V, Dec, Ra
# due to the maximum length of the data array and blocks,
# there zeros on the end of some rows and columns
data_arr = numpy.ndarray(buffer=data_buf, dtype='d', shape=(int(7*max_block_num),max_block_length))
# number of galaxy in each block
label_arr = numpy.ndarray(buffer=label_buf, dtype='i', shape=(max_block_num,max_block_length))

# loop the areas
for area in areas:
    # open the hdf5 file and assign to the array in memory
    set_path = "/grid/%s"%area
    h5_file = h5py.File(data_h5_path,"r")
    # the attribute, "grid_info", of each area
    # is a (1,2) array, [grid_rows, grid_cols].
    grid_rows, grid_cols = h5_file[set_path].attrs["grid_info"][0]
    grid_num = int(grid_rows*grid_cols)
    block_list = [i for i in range(grid_num)]
    sub_block_list = tool_box.allot(block_list, cpus)[rank]

    # the (n,4) array of the n blocks in each grid
    block_boundy = h5_file[set_path + "/boundy"].value
    block_boundx = h5_file[set_path + "/boundx"].value
    # assign to the array in memory

    if rank == 0:
        for i in range(grid_num):
            irow, icol = divmod(i, grid_cols)
            # the data in each block is a (7, n) array
            # G1, G2, N, U, V, Ra, Dec
            sub_data = h5_file[set_path+"/%d/%d"%(irow, icol)].value
            y,x = h5_file[set_path+"/%d/%d"%(irow, icol)].attrs["shape"][0]
            data_arr[int(i*7):int((1+1)*7),0:x] = sub_data
            label_arr[i,0] = int(y)
    comm.Barrier()
    h5_file.close()

    # loop the radius scales
    for ir in range(scale_num):
        radius_s, radius_e = radius_scales[ir], radius_scales[ir+1]
        # loop the mission blocks
        for ib in sub_block_list:
            # from 1-d block labels to 2-d block coordinates
            iy, ix = divmod(ib, grid_cols)
            # loop the galaxies in the block (iy, ix)
            # the shape of the data in block (iy, ix) should be stored in the attributes
            for ig in range(label_arr[ib,0]):
                dec, ra = data_arr[int(ib*7+5),ig],data_arr[int(ib*7+6),ig]
                target_blocks = tool_box.find_block(block_scale, radius_s, radius_e, iy, ix, dec, ra,
                                                    grid_rows, grid_cols, block_boundy, block_boundx)
                for itb in target_blocks:


            block_data_h5_path = block_scale_h5_path + "%d/%d/block"%(iy,ix)
            gal_num = h5_data[block_data_h5_path].attrs['shape'][0]
            block_corner_y, block_corner_x = h5_data[block_data_h5_path].attrs['origin']

