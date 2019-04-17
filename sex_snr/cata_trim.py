import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
from subprocess import Popen
import os
import time
from mpi4py import MPI
import tool_box
import numpy
from sys import argv
import h5py



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

filter_name = ["sex2_4", "sex2_2", "sex2_1.5", "sex3_4", "sex3_2", "sex3_1.5","sex4_4", "sex4_2", "sex4_1.5"]

source = argv[1]

ini_path = "%s/work/envs/envs.dat"%my_home
total_path = tool_box.config(ini_path, ['get'], [['selection_bias', "%s_path"%source, '1']])[0]


snr_idx = 0
flux_auto_idx = 1
flux_err_idx = 2
mag_auto_idx = 3
area_idx = 4
x_idx = 5
y_idx = 6


for f_nm in filter_name:

    t1 = time.time()

    data_path = total_path + "result/data/%s/sex_%d.npz"%(f_nm, rank)
    data = numpy.load(data_path)['arr_0']

    final_data = numpy.zeros((data.shape[0], 4))
    # print(final_data.shape, data[:, mag_auto_idx].shape, data.shape, len(data))

    # magnitude
    final_data[:, 0] = -data[:, mag_auto_idx]

    # snr
    final_data[:, 1] = data[:, snr_idx]

    mask = numpy.ones((len(data), 1), dtype=numpy.intc)
    idx = data[:,snr_idx] <= 0
    mask[idx] = 0

    # snr_auto
    data[:,flux_err_idx][idx] = 1
    final_data[:, 2] = data[:,flux_auto_idx]/data[:,flux_err_idx]

    # area
    final_data[:, 3] = data[:,area_idx]

    h5path = total_path + "result/data/%s/sex_%d.hdf5"%(f_nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = final_data
    h5f.close()

    h5path = total_path + "result/data/%s/mask_%d.hdf5"%(f_nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = mask
    h5f.close()

    t2 = time.time()
    print(t2-t1, mask.sum(), f_nm)
