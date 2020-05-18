import h5py
from sys import argv


data_path = argv[1]
shear_num = int(argv[2])
noisy_nm = argv[3]
noise_free_nm = argv[4]
noise_residual_nm = argv[5]
cross_term_nm = argv[6]

for i in range(shear_num):
    h5f = h5py.File(data_path + "/data_%s_%d.hdf5"%(noisy_nm,i),"r")
    data_n = h5f["/data"][()]
    h5f.close()

    h5f = h5py.File(data_path + "/data_%s_%d.hdf5"%(noise_free_nm,i),"r")
    data_nf = h5f["/data"][()]
    h5f.close()

    h5f = h5py.File(data_path + "/data_%s_%d.hdf5"%(cross_term_nm,i),"r")
    data_ct = h5f["/data"][()]
    h5f.close()

    h5f = h5py.File(data_path +"/data_%s_%d.hdf5"%(noise_residual_nm,i),"r")
    data_nr = h5f["/data"][()]
    h5f.close()

    diff = data_n - data_nf - data_ct - data_nr
    print(diff.mean(), diff.min(), diff.max())