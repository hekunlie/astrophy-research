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


# separate the data into pieces for conve

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

# the directory
total_path = argv[1]

filter_name = ["sex2_4", "sex2_2", "sex2_1.5", "sex4_4", "sex4_2", "sex4_1.5"]
# "sex3_4", "sex3_2", "sex3_1.5",


snr_idx = 0
flux_auto_idx = 1
flux_err_idx = 2
mag_auto_idx = 3
area_idx = 4
x_idx = 5
y_idx = 6


fourier_path = "%s/result/data/data_1.5sig/data_%d.hdf5"%(total_path,rank)
h5f = h5py.File(fourier_path,"r")
fourier_data = h5f["/data"][()]
h5f.close()

pk0 = fourier_data[:,4]
pk0_fit = fourier_data[:,5]
pk0_fit_max = fourier_data[:,6]

for nm in filter_name:

    t1 = time.time()

    data_path = total_path + "/result/data/%s/sex_%d.npz"%(nm, rank)
    data = numpy.load(data_path)['arr_0']

    final_data = numpy.zeros((data.shape[0], 4))
    # print(final_data.shape, data[:, mag_auto_idx].shape, data.shape, len(data))

    mask = numpy.ones((len(data), ), dtype=numpy.intc)
    idx = data[:,snr_idx] <= 0
    mask[idx] = 0

    # true magnitude
    h5path = total_path + "/parameters/para_%d.hdf5"%rank
    h5f = h5py.File(h5path,"r")
    mag_true = h5f["/mag"].value
    mag_true.shape = (mag_true.shape[0],)
    h5f.close()
    h5path = total_path + "/result/data/%s/mag_true_%d.hdf5"%(nm,rank)
    h5f = h5py.File(h5path,"w")
    h5f["/data"] = -mag_true
    h5f.close()

    # magnitude
    h5path = total_path + "/result/data/%s/mag_auto_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = -data[:, mag_auto_idx]
    h5f.close()

    # snr
    h5path = total_path + "/result/data/%s/snr_sex_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = data[:, snr_idx]
    h5f.close()

    # snr_auto
    data[:,flux_err_idx][idx] = 1
    h5path = total_path + "/result/data/%s/snr_auto_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = data[:,flux_auto_idx]/data[:,flux_err_idx]
    h5f.close()

    # area
    h5path = total_path + "/result/data/%s/area_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = data[:,area_idx]
    h5f.close()

    # Pk0
    h5path = total_path + "/result/data/%s/flux2_ex1_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = pk0
    h5f.close()

    # Pk0_fit
    h5path = total_path + "/result/data/%s/flux2_ex2_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = pk0_fit
    h5f.close()

    # max(Pk0,Pk0_fit)
    h5path = total_path + "/result/data/%s/flux2_ex3_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = pk0_fit_max
    h5f.close()

    # mask
    h5path = total_path + "/result/data/%s/mask_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = mask
    h5f.close()

    t2 = time.time()
    print(t2-t1, mask.sum(), nm)
