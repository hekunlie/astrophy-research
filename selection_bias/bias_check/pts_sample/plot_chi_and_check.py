import h5py
from mpi4py import MPI
import numpy
from sys import argv, path
path.append("D:/Github/astrophy-research/mylib")
path.append("/home/hklee/work/mylib")
import os
from plot_tool import Image_Plot


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

cmd = argv[1]

parent_path = "/mnt/perc/hklee/bias_check/data_from_pi/pow_noise_test"
src_path = parent_path + '/data_4'

dst = ["noise_free", "cross_term", "noise_residual", "noise_residual_cross_term",
       "pgal_cross_term", "pgal_noise_residual", "noisy"]

src = [["noise_free_epsf"], ["cross_term_epsf"], ["noise_diff_epsf"],
       ["cross_term_epsf", "noise_diff_epsf"], ["noise_free_epsf","cross_term_epsf"],
       ["noise_free_epsf","noise_diff_epsf"], ["noisy_cpp_epsf"]]


if cmd == "plot":
    img = Image_Plot()
    img.subplots(2,5)
    for i in range(2):
        for j in range(5):
            m = i*5 + j

            h5f = h5py.File(parent_path + "/%s/chi_%d_epsf.hdf5" % (dst[rank], m), "r")
            chisq1 = h5f["/chisq1"][()]
            chisq2 = h5f["/chisq2"][()]
            g = h5f["/g"][()]
            h5f.close()
            img.axs[i][j].plot(g, chisq1, label="g1")
            img.axs[i][j].plot(g, chisq2, label="g2")
            img.axs[i][j].legend()
            img.set_label(i,j,0,"$\chi^2$")
            img.set_label(i,j,0,"$g$")

    pic_nm = parent_path + "/%s/chisq_%s.png"%(dst[rank],dst[rank])
    img.save_img(pic_nm)


if cmd == "check":

    diff = numpy.zeros((1, 10))
    print(dst[rank], src[rank])
    for i in range(10):
        h5f = h5py.File(parent_path + "/%s/data_%d_epsf.hdf5"%(dst[rank],i),"r")
        data = h5f["/data"][()]
        h5f.close()

        for ss in src[rank]:
            h5f = h5py.File(src_path + "/data_%d_%s.hdf5" % (i, ss), "r")
            temp = h5f["/data"][()]
            h5f.close()
            data -= temp
        diff[0,i] = numpy.max(numpy.abs(data))

    print(dst[rank], src[rank],diff)

