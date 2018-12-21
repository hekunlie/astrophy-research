import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/' % my_home)
import numpy
from Fourier_Quad import Fourier_Quad
import time
from mpi4py import MPI
import h5py
import tool_box
import matplotlib.pyplot as plt

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

source = "debug1"
envs_path = "%s/work/envs/envs.dat" % my_home
get_contents = [['selection_bias', "%s_path_result" % source, '1'],['selection_bias', "%s_path_para" % source, '1']]
path_items = tool_box.config(envs_path, ['get',"get"], get_contents)
result_path,para_path = path_items

shear_cata = para_path + "shear.npz"
shear = numpy.load(shear_cata)
g1 = shear["arr_0"][:,0]
g2 = shear["arr_1"][:,0]

seed = rank*12443+1654
rng = numpy.random.RandomState(seed)
shear_num = 10
gal_num = [10, 30, 50, 200, 500]
fq = Fourier_Quad(60, 152356)

if rank == 0:
    pool = []


for num in gal_num:
    t1 = time.time()
    md = numpy.zeros((4, shear_num))
    for i in range(shear_num):

        data_path = result_path + "data/data_1.5sig/data_%d.hdf5"%i
        h5f = h5py.File(data_path,"r")
        detect = h5f["/data"].value[:,1]
        idx = detect > 0
        h5f.close()

        data_path = result_path + "data/data_%d.hdf5"%i
        h5f = h5py.File(data_path,"r")
        data = h5f["/data"].value
        h5f.close()
        MG1 = data[:, 2][idx]
        MG2 = data[:, 3][idx]
        MN = data[:, 4][idx]
        MU = data[:, 5][idx]
        MV = data[:, 6][idx]
        # be careful that the "MU" defined in FRESH is different from that in ours
        # MN + (-) MU for our definition (g1(2)) of MU and MV which is the same as
        # those in the paper Zhang et al. 2017 ApJ, 834:8
        DE1 = MN + MU
        DE2 = MN - MU
        tag = [i for i in range(len(MG1))]
        if rank == 0:
            print(len(MG1))
        ch = rng.choice(tag, num * 10000)
        g1_h, g1_sig = fq.fmin_g_new(MG1[ch], DE1[ch], bin_num=8)
        g2_h, g2_sig = fq.fmin_g_new(MG2[ch], DE2[ch], bin_num=8)
        md[0,i],md[1,i],md[2,i],md[3,i] = g1_h, g1_sig, g2_h, g2_sig

    e1mc = tool_box.data_fit(g1, md[0], md[1])
    e2mc = tool_box.data_fit(g2, md[2], md[3])
    plt.figure(figsize=(16,8))
    plt.subplot(121)
    lb = "m: %.5f(%.5f), \nc: %.6f(%.6f)"%(e1mc[0]-1,e1mc[1], e1mc[2], e1mc[3])
    plt.errorbar(g1, md[0], md[1],fmt="none",label=lb)
    plt.plot([-0.04, 0.04],[-0.04,0.04])
    plt.legend()
    plt.subplot(122)
    lb = "m: %.5f(%.5f), \nc: %.6f(%.6f)" % (e2mc[0] - 1, e2mc[1], e2mc[2], e2mc[3])
    plt.errorbar(g2, md[2], md[3],fmt="none",label=lb)
    plt.plot([-0.04, 0.04],[-0.04,0.04])
    plt.legend()
    plt.savefig(result_path+"pic/%d_%d.png"%(num,rank))
    plt.close()

    comm.Barrier()
    t2 = time.time()
    e1mcs = comm.gather(e1mc,root=0)
    e2mcs = comm.gather(e2mc,root=0)
    if rank == 0:
        e1mcs = numpy.array(e1mcs)
        e2mcs = numpy.array(e2mcs)
        numpy.savez(result_path+"pic/%d.npz"%num, e1mcs, e2mcs)
        pool.append((e1mcs, e2mcs))

    comm.Barrier()

if rank == 0:

    # for num in gal_num:
    #     a = numpy.load(result_path + "pic/%d.npz" % num)
    #     pool.append((a["arr_0"], a["arr_1"]))

    y = numpy.ones((cpus,))
    a = 1
    b = 12
    cols = ['C0', 'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9']
    for i in range(cpus):
        a += 0.6
        y[i] = a

    plt.figure(figsize=(18,18))
    plt.subplot(221)
    for i in range(len(gal_num)):
        mc1 = pool[i][0]
        lb = "$m_1$, num: $%d * 10^4$"%gal_num[i]
        plt.errorbar(x=mc1[:,0] - 1, y=y+i*b,xerr=mc1[:,1],fmt="none",label=lb,capsize=3, c=cols[i])
        plt.scatter(x=mc1[:, 0] - 1, y=y + i * b, color=cols[i])
    y1,y2 = plt.ylim()
    plt.plot([-0.01, -0.01], [y1, y2], c="grey", linestyle="--")
    plt.plot([0.01, 0.01], [y1, y2], c="grey", linestyle="--")
    plt.plot([0, 0],[y1,y2],c="grey",linestyle="--")
    plt.ylim(y1,y2)
    plt.legend()
    plt.subplot(222)
    for i in range(len(gal_num)):
        mc1 = pool[i][0]
        lb = "$c_1$, num: $%d * 10^4$"%gal_num[i]
        plt.errorbar(x=mc1[:,2],y=y+i*b,xerr=mc1[:,3],fmt="none",label=lb,capsize=3,c=cols[i])
        plt.scatter(x=mc1[:,2],y=y+i*b,color=cols[i])
    y1,y2 = plt.ylim()
    plt.plot([0,0],[y1,y2],c="grey",linestyle="--")
    plt.plot([-0.0001, -0.0001], [y1, y2], c="grey", linestyle="--")
    plt.plot([0.0001, 0.0001], [y1, y2], c="grey", linestyle="--")
    plt.ylim(y1,y2)
    plt.legend()
    plt.savefig(result_path + "pic/total.png")

    plt.subplot(223)
    for i in range(len(gal_num)):
        mc2 = pool[i][1]
        lb = "$m_2$, num: $%d * 10^4$"%gal_num[i]
        plt.errorbar(x=mc2[:,0]-1,y=y+i*b,xerr=mc2[:,1],fmt="none",label=lb,capsize=3, c=cols[i])
        plt.scatter(x=mc2[:, 0]-1, y=y + i * b, color=cols[i])
    y1,y2 = plt.ylim()
    plt.plot([0,0],[y1,y2],c="grey",linestyle="--")
    plt.plot([-0.01, -0.01], [y1, y2], c="grey", linestyle="--")
    plt.plot([0.01, 0.01], [y1, y2], c="grey", linestyle="--")
    plt.ylim(y1,y2)
    plt.legend()
    plt.subplot(224)
    for i in range(len(gal_num)):
        mc2 = pool[i][1]
        lb = "$c_2$, num: $%d * 10^4$"%gal_num[i]
        plt.errorbar(x=mc2[:,2],y=y+i*b,xerr=mc2[:,3],fmt="none",label=lb,capsize=3, c=cols[i])
        plt.scatter(x=mc2[:, 2], y=y + i * b, color=cols[i])
    y1,y2 = plt.ylim()
    plt.plot([0,0],[y1,y2],c="grey",linestyle="--")
    plt.plot([-0.0001, -0.0001], [y1, y2], c="grey", linestyle="--")
    plt.plot([0.0001, 0.0001], [y1, y2], c="grey", linestyle="--")
    plt.ylim(y1,y2)
    plt.legend()
    plt.savefig(result_path+"pic/total.png")