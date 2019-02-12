import numpy
import h5py
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

mu = [0, 0]
data_len = 5000000
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
data_path = "/mnt/ddnfs/data_users/hkli/correlation/simu/mgauss_%d.hdf5"%rank
h5f = h5py.File(data_path,"w")
h5f["/data"] = data
h5f.close()