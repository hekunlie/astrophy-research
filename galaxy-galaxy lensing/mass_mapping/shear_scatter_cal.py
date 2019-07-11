import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import numpy
import h5py
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import tool_box
import matplotlib.pyplot as plt
from mpi4py import MPI
import time


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cmd = argv[1]

parent_path = "/mnt/perc/hklee/CFHT/multi_shear/result/"
h5f = h5py.File("/mnt/perc/hklee/CFHT/multi_shear/data/shear.hdf5", "r")
g1 = h5f["/g1"].value
g2 = h5f["/g2"].value
h5f.close()

shear_num = g1.shape[0]*g1.shape[1]

g1 = g1.reshape((shear_num,))
g2 = g2.reshape((shear_num,))

gal_each_shear = 500000
num = [1, 5, 10, 20, int(gal_each_shear/10000)]
ch_num = [10000, 50000, 100000, 200000, gal_each_shear]


if cmd == "cal":

    t1 = time.time()

    itemsize = MPI.DOUBLE.Get_size()
    if rank == 0:
        nbytes_data = shear_num*gal_each_shear*5*itemsize
        nbytes = (shear_num*2+2)*itemsize*5
    else:
        nbytes_data = 0
        nbytes = 0

    win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
    win2 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
    win3 = MPI.Win.Allocate_shared(nbytes_data, itemsize, comm=comm)

    buf1, itemsize = win1.Shared_query(0)
    buf2, itemsize = win2.Shared_query(0)
    buf3, itemsize = win3.Shared_query(0)

    result_1 = numpy.ndarray(buffer=buf1, dtype='d', shape=(len(ch_num), shear_num*2 + 2)) # g1
    result_2 = numpy.ndarray(buffer=buf2, dtype='d', shape=(len(ch_num), shear_num*2 + 2)) # g2
    total_data_buf = numpy.ndarray(buffer=buf3, dtype='d', shape=(shear_num*gal_each_shear, 5)) # g2

    my_task = tool_box.allot([i for i in range(shear_num)], cpus)[rank]
    print("Rank: %d. "%rank, my_task)
    fq = Fourier_Quad(12,11)
    for ig in my_task:
        data_path = parent_path + "shear_%d.hdf5"%ig
        h5f = h5py.File(data_path,"r")
        data = h5f["/data"].value
        h5f.close()

        mg1 = data[:,0]
        mg2 = data[:,1]
        mnu1 = data[:,2] + data[:,3]
        mnu2 = data[:,2] - data[:,3]
        total_data_buf[ig*gal_each_shear:(ig+1)*gal_each_shear] = data
        for i in range(len(ch_num)):
            gh1, gh1_sig = fq.fmin_g_new(g=mg1[:ch_num[i]], nu=mnu1[:ch_num[i]], bin_num=8, scale=100)[:2]
            gh2, gh2_sig = fq.fmin_g_new(g=mg2[:ch_num[i]], nu=mnu2[:ch_num[i]], bin_num=8, scale=100)[:2]

            result_1[i, ig] = gh1
            result_1[i, ig + shear_num] = gh1_sig
            result_2[i, ig] = gh2
            result_2[i, ig + shear_num] = gh2_sig

    comm.Barrier()
    t2 = time.time()
    if rank == 0:

        for j in range(shear_num):
            h5f = h5py.File(parent_path + "shear_%d.hdf5" % j, "r")
            temp = h5f["/data"].value
            h5f.close()
            if j == 0:
                data = temp
            else:
                data = numpy.row_stack((data, temp))

        img = Image_Plot()
        img.subplots(2, len(ch_num))
        img.set_style()

        for i in range(len(ch_num)):

            idx = [n + m*gal_each_shear for m in range(shear_num) for n in range(ch_num[i])]

            mg1 = data[:, 0][idx]
            mg2 = data[:, 1][idx]
            mnu1 = data[:, 2][idx] + data[:, 3][idx]
            mnu2 = data[:, 2][idx] - data[:, 3][idx]

            gh1, gh1_sig = fq.fmin_g_new(g=mg1, nu=mnu1, bin_num=8, scale=100,fig_ax=img.axs[0][i])[:2]
            gh2, gh2_sig = fq.fmin_g_new(g=mg2, nu=mnu2, bin_num=8, scale=100,fig_ax=img.axs[1][i])[:2]

            result_1[i, shear_num*2] = gh1
            result_1[i, shear_num*2 + 1] = gh1_sig
            result_2[i, shear_num*2] = gh2
            result_2[i, shear_num*2 + 1] = gh2_sig
            print(g1.mean(), gh1, gh1_sig)
            print(g2.mean(), gh2, gh2_sig)
        img.save_img(parent_path + "pic/check_total.png")
        img.close_img()

        h5f = h5py.File(parent_path + "final_result.hdf5","w")
        h5f["/g1"] = result_1
        h5f["/g2"] = result_2
        h5f.close()

        t3 = time.time()
        print("Done: %.2f sec, %.2f sec"%(t2-t1, t3-t2))


if cmd == "plot":

    h5f = h5py.File(parent_path + "final_result.hdf5", "r")
    result_1 = h5f["/g1"].value
    result_2 = h5f["/g2"].value
    h5f.close()
    print(g1.mean())
    print(g2.mean())
    # check
    img = Image_Plot(fig_x=12, fig_y=9)
    img.subplots(1,2)
    img.set_style()
    for j in range(len(ch_num)):

        img.axs[0][0].errorbar(g1, result_1[j, :shear_num], result_1[j, shear_num:2*shear_num], capsize=5, marker=">",
                               fmt="none", c="C%d"%j, label="$NUM: %d \\times 10^4$"%num[j])
        img.axs[0][0].scatter(g1, result_1[j, :shear_num], marker=">", c="C%d"%j)
        img.axs[0][0].errorbar(g1.mean(), result_1[j, 2*shear_num], result_1[j, 2*shear_num+1], capsize=5, marker="s",
                               fmt="none", c="C%d"%j)
        img.axs[0][0].scatter(g1.mean(), result_1[j, 2*shear_num],marker="s",c="C%d"%j)
        img.axs[0][0].plot(g1, g1, c="black")

        img.axs[0][1].errorbar(g2, result_2[j, :shear_num], result_2[j, shear_num:2*shear_num], capsize=5, marker=">",
                               fmt="none", c="C%d"%j, label="$NUM: %d \\times 10^4$"%num[j])
        img.axs[0][1].scatter(g2, result_2[j, :shear_num],marker=">", c="C%d"%j)
        img.axs[0][1].errorbar(g2.mean(), result_2[j, 2*shear_num], result_2[j, 2*shear_num+1], capsize=5, marker="s",
                               fmt="none", c="C%d"%j)
        img.axs[0][1].scatter(g2.mean(), result_2[j, 2*shear_num],marker="s",c="C%d"%j)
        img.axs[0][1].plot(g2, g2, c="black")

    img.axs[0][0].legend()
    img.axs[0][1].legend()
    img.set_label(0,0,0,"$g1_{m}$")
    img.set_label(0,0,1,"$g1_{true}$")
    img.set_label(0,1,0,"$g2_{m}$")
    img.set_label(0,1,1,"$g2_{true}$")
    img.save_img(parent_path + "final_result_check.png")
    img.close_img()

    # plot the error bar vs gal_num
    img = Image_Plot()
    img.subplots(1, 1)

    sigma = result_1[:,shear_num:2*shear_num].sum(axis=1)/shear_num

    img.axs[0][0].plot(ch_num, sigma, c="C1")
    img.axs[0][0].plot(ch_num, result_1[:, shear_num*2+1], c="C2")
    img.axs[0][0].plot(ch_num, result_1[:, shear_num * 2 + 1]*numpy.sqrt(shear_num), c="C4")
    img.axs[0][0].scatter(ch_num, sigma, c="C1", label="$Mean \ \sigma$")
    img.axs[0][0].scatter(ch_num, result_1[:,shear_num*2+1], c="C2", label="$\sigma\  of\  all$")
    img.axs[0][0].scatter(ch_num, result_1[:,shear_num*2+1]*numpy.sqrt(shear_num), c="C4", label="$real\ \sigma\  of\  all$")
    img.axs[0][0].legend()
    img.axs[0][0].set_xscale("log")
    img.set_label(0,0,0,"$\sigma$ of $g_1$")
    img.set_label(0,0,1,"$Num$")

    img.save_img(parent_path + "final_result_g1.png")
    img.close_img()

    img = Image_Plot()
    img.subplots(1, 1)

    sigma = result_2[:,shear_num:2*shear_num].sum(axis=1)/shear_num

    img.axs[0][0].plot(ch_num, sigma, c="C1")
    img.axs[0][0].plot(ch_num, result_2[:,shear_num*2+1], c="C2")
    img.axs[0][0].plot(ch_num, result_2[:, shear_num * 2 + 1]*numpy.sqrt(shear_num), c="C4")
    img.axs[0][0].scatter(ch_num, sigma, c="C1", label="$Mean \ \sigma$")
    img.axs[0][0].scatter(ch_num, result_2[:,shear_num*2+1], c="C2", label="$\sigma\  of\  all$")
    img.axs[0][0].scatter(ch_num, result_2[:,shear_num*2+1]*numpy.sqrt(shear_num), c="C4", label="$real\  \sigma\  of\  all$")
    img.axs[0][0].legend()
    img.axs[0][0].set_xscale("log")
    img.set_label(0,0,0,"$\sigma$ of $g_2$")
    img.set_label(0,0,1,"$Num$")

    img.save_img(parent_path + "final_result_g2.png")
    img.close_img()