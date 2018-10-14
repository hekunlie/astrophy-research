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


# the data and directories have not been prepared !!!

areas = ['w1', 'w2', 'w3', 'w4']
scale_num = 15
radius_scales = [0]
for i in range(scale_num):
    radius_scales.append(10**(0.15 * i))

# loop the areas in each h5-file
for area in areas:
    # open the hdf5 file which contains all the data in the mesh grid
    data_h5_path = data_path + ".../"
    h5_data = h5py.File(data_h5_path)

    # loop the radius scales
    for ir in range(scale_num):
        radius_s, radius_e = radius_scales[ir], radius_scales[ir+1]

        # find the block scale needed, the number of blocks should be stored in the attributes
        block_scale = find_block_scale()!!!

        # a tuple of shape of the block, (y,x)
        block_scale_h5_path = "/%d/"%block_scale
        block_ny, block_nx = h5_data[block_scale_h5_path].attrs['shape'] !!!
        corner_y, corner_x = h5_data[block_scale_h5_path].attrs['origin']!!!

        block_labels = [i for i in range(block_nx*block_ny)]
        mission_blocks = tool_box.allot(block_labels, cpus)[rank]

        # loop the mission blocks
        for ig in mission_blocks:
            # from 1-d block labels to 2-d block coordinates
            iy, ix = divmod(ig, block_nx)
            # loop the galaxies in the block (iy, ix)
            # the shape of the data in block (iy, ix) should be stored in the attributes
            block_data_h5_path = block_scale_h5_path + "%d/%d/block"%(iy,ix)
            gal_num = h5_data[block_data_h5_path].attrs['shape'][0]
            block_corner_y, block_corner_x = h5_data[block_data_h5_path].attrs['origin']

