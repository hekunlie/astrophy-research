import matplotlib
matplotlib.use("Agg")
import numpy
import matplotlib.pyplot as plt
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import time
from Fourier_Quad import Fourier_Quad
from sys import argv
from mpi4py import MPI
import tool_box
import h5py

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cut = argv[1]

t1 = time.clock()

shear = numpy.load("/home/hklee/work/selection_bias/parameters/shear.npz")["arr_0"]
fg1 = shear[0]
fg2 = shear[1]

num = 10000000
sp = (10000000, 27)
if rank == 0:
    arr = numpy.load("/lmc/selection_bias/result_kw_new/data/data_cache.npz")["arr_0"]
    #tag = arr[:, 4]
    data = arr[0: num, :]
    for i in range(1, 10):
        data_s = arr[i*num: (i+1)*num, :]
        comm.Send(data_s, dest=i, tag=i)

else:
    data = numpy.empty(sp, dtype=numpy.float64)
    comm.Recv(data, source=0, tag=rank)


fq = Fourier_Quad(52, 152356)
nsig = 380.64
area = data[:, 16]

flux = data[:, 17]/nsig
fcut = [0, 80, 100, 115, 125, 135, 145, 155, 200, 300, 450]

peak = data[:, 18]/nsig
pcut = [0, 5, 6, 7, 7.5, 8, 9, 10, 12, 17, 22]

snr = data[:, 19]
scut = [0, 18, 21, 23, 25, 27, 29, 35, 40, 50, 65]

fsnr = data[:, 20]
fscut = [0, 2.8, 3.3, 3.8, 4.1, 4.5, 4.9, 5.3, 7, 10, 15]

osnr = data[:, 21]
ocut = [0, 18, 21, 23, 25, 27, 29, 35, 40, 50, 65]

N = data[:10000000, 10]
n = numpy.sort(N)
ncut = [n[int(0.1*len(n)*i)] for i in range(1, 8)]

select = {"N": (N, ncut), "flux": (flux, fcut), "peak": (peak, pcut),"fsnr": (fsnr, fscut), "snr": (snr, scut), "osnr": (osnr, ocut)}

res_arr = numpy.zeros((6, len(select[cut][1])))

mg1 = data[:, 3]
mg2 = data[:, 8]
mn = data[:, 10]
mu = data[:, 11]
mv = data[:, 12]

for tag, cut_s in enumerate(select[cut][1]):
    idx = select[cut][0] > cut_s
    num = len(mg1[idx])

    g1_h, g1_sig = fq.fmin_g(mg1[idx], mn[idx], mu[idx], mode=1, bin_num=8)
    g2_h, g2_sig = fq.fmin_g(mg2[idx], mn[idx], mu[idx], mode=2, bin_num=8)

    res_arr[:, tag] = numpy.array([g1_h, g1_sig, num, g2_h, g2_sig, num])

if rank > 0:
    comm.Send(res_arr, dest=0, tag=rank)

else:
    for procs in range(1, 10):
        recvs = numpy.empty((6, len(select[cut][1])), dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        res_arr = numpy.column_stack((res_arr, recvs))

    mc1 = []
    mc2 = []
    for tag, cut_s in enumerate(select[cut][1]):
        arr = res_arr[:, tag]
        for i in range(1, 10):
            arr = numpy.column_stack((arr, res_arr[:, i*len(select[cut][1])+tag]))

        e1mc = tool_box.data_fit(fg1, arr[0], arr[1])
        mc1.append(e1mc)
        e2mc = tool_box.data_fit(fg2, arr[3], arr[4])
        mc2.append(e2mc)

        mc = numpy.array([e1mc, e2mc])
        data_path = "/home/hklee/work/result/cuts/sym/" + cut + "/" + str(cut_s)+".npz"
        numpy.savez(data_path, arr, mc)
        pic_path = "/home/hklee/work/result/cuts/sym/" + cut + "/" + str(cut_s)+".eps"
        tool_box.mcplot(fg1, arr[0:3,:], fg2, arr[3:6,:], e1mc, e2mc, str(cut_s), 'max', pic_path)
        pic_path = "/home/hklee/work/result/cuts/sym/" + cut + "/" + str(cut_s) + ".png"
        tool_box.mcplot(fg1, arr[0:3, :], fg2, arr[3:6, :], e1mc, e2mc, str(cut_s), 'max', pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = "/home/hklee/work/result/cuts/sym/" + cut + "/total.npz"
    numpy.savez(mc_path, mc1, mc2)

    x1 = select[cut][1][0]
    x2 = select[cut][1][-1]

    fig = plt.figure(figsize=(10, 10))

    plt.subplot(221)
    plt.errorbar(select[cut][1], mc1[0] - 1, mc1[1], c='coral', capsize=2)
    #plt.plot([0.05 * (x1 - x2), 1.05 * (x2 - x1)], [0, 0], c='grey')
    #plt.xlim(0.05 * (x1 - x2), 1.05 * (x2 - x1))
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((1, 2))
    plt.xlabel("Cutoff")
    plt.semilogx()
    plt.ylabel("m1")

    plt.subplot(222)
    plt.errorbar(select[cut][1], mc1[2], mc1[3], c='coral', capsize=2)
   # plt.plot([0.05 * (x1 - x2), 1.05 * (x2 - x1)], [0, 0], c='grey')
    #plt.xlim(0.05 * (x1 - x2), 1.05 * (x2 - x1))
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((1, 2))
    plt.xlabel("Cutoff")
    plt.semilogx()
    plt.ylabel("c1")

    plt.subplot(223)
    plt.errorbar(select[cut][1], mc2[0] - 1, mc2[1], c='coral', capsize=2)
    #plt.plot([0.05 * (x1 - x2), 1.05 * (x2 - x1)], [0, 0], c='grey')
    #plt.xlim(0.05 * (x1 - x2), 1.05 * (x2 - x1))
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((1, 2))
    plt.xlabel("Cutoff")
    plt.semilogx()
    plt.ylabel("m2")

    plt.subplot(224)
    plt.errorbar(select[cut][1], mc2[2], mc2[3], c='coral', capsize=2)
    #plt.plot([0.05 * (x1 - x2), 1.05 * (x2 - x1)], [0, 0], c='grey')
    #plt.xlim(0.05 * (x1 - x2), 1.05 * (x2 - x1))
    ax = plt.gca()
    ax.yaxis.get_major_formatter().set_powerlimits((1, 2))
    plt.xlabel("Cutoff")
    plt.semilogx()
    plt.ylabel("c2")

    namep = "/home/hklee/work/result/cuts/sym/" + cut + "/total.eps"
    plt.savefig(namep)
    plt.show()
    plt.close()

t2 = time.clock()
if rank == 0:
    print(t2-t1)

