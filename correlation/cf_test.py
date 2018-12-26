import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import time
from mpi4py import MPI
from Fourier_Quad import Fourier_Quad
import tool_box
import numpy
import h5py
from sys import argv


sect, source = argv[1], argv[2]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['%s'%sect, "%s_path" % source, '1'], ['%s'%sect, "%s_path_result" % source, '1'],
                ['%s'%sect, "%s_path_para" % source, '1'], ['%s'%sect, "%s_path_log" % source, '1']]
path_items = tool_box.config(envs_path, ['get', 'get', 'get', 'get'], get_contents)
total_path, result_path, para_path, log_path = path_items

logger = tool_box.get_logger(log_path + "%d_logs.dat" % rank)

# the parameters
para_contents = [["para", "total_num", 1], ["para", "stamp_size", 1], ["para", "stamp_col", 1], ["para", "shear_num", 1],
                 ["para", "noise_sig", 1], ["para", "pixel_scale", 1]]
para_items = tool_box.config(para_path+"para.ini", ['get', 'get', 'get', 'get', 'get', 'get'], para_contents)

total_chips_num = int(para_items[0])
stamp_size = int(para_items[1])
stamp_col = int(para_items[2])
shear_num = int(para_items[3])
noise_sig = int(para_items[4])
pixel_scale = float(para_items[5])
stamp_num = 10000

fq = Fourier_Quad(64, 123)
cov11, cov12 = 0.0006, -0.0003
if rank < int(shear_num/2):
    h5path = result_path + "/data/data_%d.hdf5"%rank
    h5f = h5py.File(h5path,"r")
    data_1 = h5f["/data"].value
    h5f.close()

    h5path = result_path + "/data/data_%d.hdf5"%(rank+int(shear_num/2))
    h5f = h5py.File(h5path,"r")
    data_2 = h5f["/data"].value
    h5f.close()

    mg1_1 = data_1[:,2]
    mn1_1 = data_1[:,4] + data_1[:,5]
    mg2_1 = data_1[:, 3]
    mn2_1 = data_1[:,4] - data_1[:,5]

    mg1_2 = data_2[:,2]
    mn1_2 = data_2[:,4] + data_2[:,5]
    mg2_2 = data_2[:,3]
    mn2_2 = data_2[:,4] - data_2[:,5]

    g1_corr, g1_corr_err = fq.fmin_g2d([mg1_1, mg1_2], [mn1_1, mn1_2],8)
    g2_corr, g2_corr_err = fq.fmin_g2d([mg2_1, mg2_2], [mn2_1, mn2_2],8)
    print("RNAK: %d: chi_11: %8.6f ( %10.8f ), chi_22: %8.6f (% 10.8f )\n"
          "TRUE:    chi_11: %8.6f           , chi_22: %8.6f             "
          %(rank, g1_corr, g1_corr_err, g2_corr, g2_corr_err, cov12 + rank*0.0001,-cov12 - rank * 0.0001))