import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
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


def G_bin2d(mgs, mnus, g_distri, bins):
    r"""
    to calculate the symmetry of two sets of shear estimators
    :param mgs: a two components list, two 1-D numpy arrays, contains two sets of
                first shear estimators of Fourier quad ,
    :param mnus: a two components list, two 1-D numpy arrays, contains two sets of
                shear estimators of Fourier quad (N,U), N + U for g1, N - U for g2
    :param g_distri: float, the correlation between the two sets of shear
    :param bins: a two components list, two 1-D numpy arrays, contains two bins for
                the shear estimator
    :return: chi square
    """
    # half of the bin number
    ny, nx = int((len(bins[0]) - 1) / 2), int((len(bins[1]) - 1) / 2)
    mg1 = mgs[0] - mnus[0] * g_distri[:, 0]
    mg2 = mgs[1] - mnus[1] * g_distri[:, 1]
    num_arr = numpy.histogram2d(mg1, mg2, bins)[0]
    # | arr_1 | arr_2 |
    # | arr_3 | arr_4 |
    # chi square = 0.2*SUM[(arr_2 + arr_3 - arr_1 - arr_4)**2/(arr_1 + arr_2 + arr_3 + arr_4)]
    arr_1 = num_arr[0:ny, 0:nx][:, range(ny - 1, -1, -1)]
    arr_2 = num_arr[0:ny, nx:2 * nx]
    arr_3 = num_arr[ny:2 * ny, 0:nx][range(ny - 1, -1, -1)][:, range(nx - 1, -1, -1)]
    arr_4 = num_arr[ny:2 * ny, nx:2 * nx][range(ny - 1, -1, -1)]
    chi_sq = 0.5 * numpy.sum(((arr_2 + arr_3 - arr_1 - arr_4) ** 2) / (arr_1 + arr_2 + arr_3 + arr_4))

    return chi_sq, num_arr


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

    data_path = result_path + "data/data_%d.hdf5"%(rank+5)
    h5f = h5py.File(data_path)
    data_2 = h5f["/data"].value
    h5f.close()

    data_path = total_path + "mgauss_%d.hdf5"%rank
    h5f = h5py.File(data_path)
    cor_gs = h5f["/data"].value
    h5f.close()

    fq = Fourier_Quad(123,123)

    mg1_1 = data_1[:,2]
    mg2_1 = data_1[:,3]
    mn1_1 = data_1[:,4] + data_1[:,5]
    mn2_1 = data_1[:,4] - data_1[:,5]

    mg1_2 = data_2[:,2]
    mg2_2 = data_2[:,3]
    mn1_2 = data_2[:,4] + data_2[:,5]
    mn2_2 = data_2[:,4] - data_2[:,5]

    bins = fq.set_bin(mg1_1, 6)
    print(rank, bins)
    chi_pool = [[], []]
    gs = [j * 0.00005 - 0.001 for j in range(40)]
    for i in range(40):
        # gch11 = mg1_1 - cor_gs[:data_len, 0]*mn1_1
        # gch12 = mg1_2 - cor_gs[:data_len, 1]*mn1_2
        #
        # gch21 = mg2_1 - cor_gs[:data_len, 0]*mn2_1
        # gch22 = mg2_2 - cor_gs[:data_len, 1]*mn2_2

        chi_1, nums_1 = G_bin2d([mg1_1, mg1_2], [mn1_1, mn1_2], cor_gs[i*data_len:(i+1)*data_len], [bins, bins])
        chi_2, nums_2 = G_bin2d([mg2_1, mg2_2], [mn2_1, mn2_2], cor_gs[i*data_len:(i+1)*data_len], [bins, bins])
        # print(rank, chi_1, chi_2)
        # print(rank, nums_1)
        # print(rank, nums_2)
        chi_pool[0].append(chi_1)
        chi_pool[1].append(chi_2)
    plt.scatter(gs, chi_pool[0])
    plt.scatter(gs, chi_pool[1])
    plt.savefig("/home/hkli/work/cpp/test/pic/%d.png"%rank)
    plt.close()

