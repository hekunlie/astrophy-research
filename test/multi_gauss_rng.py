import numpy
import h5py
from mpi4py import MPI
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
from Fourier_Quad import Fourier_Quad
import tool_box

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()
sect, source = "correlation", "cor"
envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['%s'%sect, "%s_path" % source, '1'], ['%s'%sect, "%s_path_result" % source, '1']]
path_items = tool_box.config(envs_path, ['get', 'get'], get_contents)
total_path, result_path = path_items

mu = [0, 0]
data_len = 5000000
if argv[1] == "rng":
    rng = numpy.random.RandomState(rank+1000*rank+1)
    for j in range(40):
        g_corr = j * 0.00005 - 0.001
        cov = [[abs(2 * g_corr), g_corr],
               [g_corr, abs(2 * g_corr)]]
        temp = rng.multivariate_normal(mu,cov,data_len)
        print(temp.shape)
        if j == 0:
            data = temp.copy()
        else:
            data = numpy.row_stack((data, temp))
    data_path = total_path + "mgauss_%d.hdf5"%rank
    h5f = h5py.File(data_path,"w")
    h5f["/data"] = data
    h5f.close()

else:
    data_path = result_path + "data/data_%d.hdf5"%rank
    h5f = h5py.File(data_path)
    data_1 = h5f["/data"].value
    h5f.close()

    data_path = result_path + "data/data_%d.hdf5"%rank+5
    h5f = h5py.File(data_path)
    data_2 = h5f["/data"].value
    h5f.close()

    data_path = total_path + "mgauss_%d.hdf5"%rank
    h5f = h5py.File(data_path)
    cor_gs = h5f["/data"].value
    h5f.close()


    gch11 = data_1[:,2] - cor_gs[:data_len, 0]*(data_1[:,4] + data_1[:,5])
    gch12 = data_2[:,2] - cor_gs[:data_len, 1]*(data_2[:,4] + data_2[:,5])

    gch21 = data_1[:,3] - cor_gs[:data_len, 0]*(data_1[:,4] - data_1[:,5])
    gch22 = data_2[:,3] - cor_gs[:data_len, 1]*(data_2[:,4] - data_2[:,5])



