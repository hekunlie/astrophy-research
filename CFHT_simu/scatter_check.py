import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import time
from mpi4py import MPI
import h5py
import tool_box
import shutil

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

source = argv[1]

envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['selection_bias', "%s_path" % source, '1'], ['selection_bias', "%s_path_result" % source, '1'],
                ['selection_bias', "%s_path_para" % source, '1'], ['selection_bias', "%s_path_log" % source, '1']]
path_items = tool_box.config(envs_path, ['get', 'get', 'get', 'get'], get_contents)
total_path, result_path, para_path, log_path = path_items

logger = tool_box.get_logger(log_path + "%d_logs.dat" % rank)

# the parameters
para_contents = [["para","total_num",1], ["para","stamp_size",1], ["para", "stamp_col", 1], ["para","shear_num",1],
                 ["para","noise_sig",1], ["para", "pixel_scale", 1]]
para_items = tool_box.config(para_path+"para.ini", ['get', 'get', 'get', 'get', 'get', 'get'], para_contents)

total_chips_num = int(para_items[0])
stamp_size = int(para_items[1])
stamp_col = int(para_items[2])
shear_num = int(para_items[3])
noise_sig = int(para_items[4])
pixel_scale = float(para_items[5])
stamp_num = 10000