import os
from sys import argv
import numpy
import h5py
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

source_nm = argv[1]

sex_nms = ["sex2_1.5","sex2_2","sex2_4","sex4_1.5","sex4_2","sex4_4"]

# total_path = "/mnt/perc/hklee/selection_bias/paper_data/%s/result/data"%source_nm
total_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/data"%source_nm

# for nm in sex_nms:
#     for i in range(10):
#         r_path = total_path + "/%s/Rfacotr_%d.hdf5"%(nm,i)
#         n_r_path = total_path + "/%s/rfactor_%d.hdf5"%(nm,i)
#         if not os.path.exists(n_r_path):
#             os.rename(r_path, n_r_path)
#
# exit()
snr_idx = 0
flux_auto_idx = 1
flux_err_idx = 2
mag_auto_idx = 3
area_idx = 4
x_idx = 5
y_idx = 6


fourier_path = "%s/data_%d.hdf5"%rank
h5f = h5py.File(fourier_path,"r")
fourier_data = h5f["/data"].value
h5f.close()

pk0 = fourier_data[:,4]
pk0_fit = fourier_data[:,5]
pk0_fit_max = fourier_data[:,6]

for nm in sex_nms:
    data_path = total_path + "result/data/%s/sex_%d.npz"%(nm, rank)
    data = numpy.load(data_path)['arr_0']

    final_data = numpy.zeros((data.shape[0], 4))
    # print(final_data.shape, data[:, mag_auto_idx].shape, data.shape, len(data))

    mask = numpy.ones((len(data), ), dtype=numpy.intc)
    idx = data[:,snr_idx] <= 0
    mask[idx] = 0

    # magnitude
    h5path = total_path + "result/data/%s/mag_auto_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = -data[:, mag_auto_idx]
    h5f.close()

    # snr
    h5path = total_path + "result/data/%s/snr_sex_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = data[:, snr_idx]
    h5f.close()

    # snr_auto
    data[:,flux_err_idx][idx] = 1
    h5path = total_path + "result/data/%s/snr_auto_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = data[:,flux_auto_idx]/data[:,flux_err_idx]
    h5f.close()

    # area
    h5path = total_path + "result/data/%s/area_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = data[:,area_idx]
    h5f.close()

    # Pk0
    h5path = total_path + "result/data/%s/flux2_ex1_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = pk0
    h5f.close()

    # Pk0_fit
    h5path = total_path + "result/data/%s/flux2_ex2_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = pk0_fit
    h5f.close()

    # max(Pk0,Pk0_fit)
    h5path = total_path + "result/data/%s/flux2_ex3_%d.hdf5"%(nm, rank)
    h5f = h5py.File(h5path, "w")
    h5f["/data"] = pk0_fit_max
    h5f.close()