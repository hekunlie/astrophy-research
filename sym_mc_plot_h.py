import matplotlib
matplotlib.use("Agg")
import numpy
import matplotlib.pyplot as plt
from sys import path
path.append('/home/hkli/work/fourier_quad/')
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
total_path = "/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/"
shear = numpy.load("/mnt/ddnfs/data_users/hkli/selection_bias_real_dimmer/parameters/shear.npz")
fg1 = shear["arr_0"]
fg2 = shear["arr_1"]
for i in range(10):
    h5path = total_path + "result/data/data_%d_%d.hdf5"%(rank,i)
    f = h5py.File(h5path, "r")
    if i == 0:
        data = f["/data"].value
    else:
        data = numpy.row_stack((data, f["/data"].value))
    f.close()

fq = Fourier_Quad(84, 152356)
nsig = 380.64

osnr = data[:, 7]
osnrcut = [0,9,13,18, 24, 30, 36, 45, 55, 70, 90, 120, 160]

flux = data[:, 8]/nsig
fcut = [0,30,65,100, 150, 200, 250, 350, 450, 550, 800, 1000, 1400]

peak = data[:, 9]/nsig
pcut = [0, 4, 5, 6, 7, 8, 10, 13, 16, 20, 25, 30]
#
fsnr = data[:, 10]
fsnrcut = [3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 20, 25, 30]

fsnr4 = data[:, 11]
fsnr4cut = [3, 4, 5, 6, 7, 8, 9, 10,11, 13, 15, 20, 25, 30]

fsnr_c = data[:, 12]
fsnrccut = [3, 4, 5, 6, 7, 8, 9, 10,11, 13, 15, 20, 25, 30]
#
fsnr_c4 = data[:, 13]
fsnrc4cut = [3, 4, 5, 6, 7, 8, 9, 10,11, 13, 15, 20, 25, 30]

snr = data[:, 14]
snrcut = [0,9,13,18, 24, 30, 36, 45, 55, 70, 90, 120, 160]

sesnr = data[:, 15]
secut = [ 13, 15, 18, 21, 25, 30, 40, 60, 90, 120]


select = {'osnr':(osnr, osnrcut),"fsnr": (fsnr, fsnrcut), "flux": (flux, fcut), "peak": (peak, pcut), "fsnr4": (fsnr4, fsnr4cut),
          "fsnr_c": (fsnr_c, fsnrccut), "fsnr_c4": (fsnr_c4, fsnrc4cut), 'snr':(snr, snrcut)}

res_arr = numpy.zeros((6, len(select[cut][1])))

mg1 = data[:, 2]
mg2 = data[:, 3]
mn = data[:, 4]
mu = data[:, 5]
mv = data[:, 6]

for tag, cut_s in enumerate(select[cut][1]):
    idx = select[cut][0] > cut_s
    num = len(mg1[idx])

    g1_h, g1_sig = fq.fmin_g(mg1[idx], mn[idx], mu[idx], mode=1, bin_num=8)
    g2_h, g2_sig = fq.fmin_g(mg2[idx], mn[idx], mu[idx], mode=2, bin_num=8)

    res_arr[:, tag] = numpy.array([g1_h, g1_sig, num, g2_h, g2_sig, num])

if rank > 0:
    # pass
    comm.Send(res_arr, dest=0, tag=rank)

else:
    for procs in range(1, cpus):
        recvs = numpy.empty((6, len(select[cut][1])), dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        res_arr = numpy.column_stack((res_arr, recvs))

    mc1 = []
    mc2 = []
    for tag, cut_s in enumerate(select[cut][1]):
        arr = res_arr[:, tag]
        for i in range(1, cpus):
            arr = numpy.column_stack((arr, res_arr[:, i*len(select[cut][1])+tag]))

        e1mc = tool_box.data_fit(fg1, arr[0], arr[1])
        mc1.append(e1mc)
        e2mc = tool_box.data_fit(fg2, arr[3], arr[4])
        mc2.append(e2mc)

        mc = numpy.array([e1mc, e2mc])
        data_path = total_path + "result/cuts/sym/" + cut + "/" + str(cut_s)+".npz"
        numpy.savez(data_path, arr, mc)

        mc_title = ['0', '0', '0', '0']
        m_r = [[e1mc[0]-1 - 2 * e1mc[1], e1mc[0]-1 + 2 * e1mc[1]], [e2mc[0]-1 - 2 * e2mc[1], e2mc[0]-1 + 2 * e2mc[1]]]
        c_r = [[e1mc[2] - 2 * e1mc[3], e1mc[2] + 2 * e1mc[3]], [e2mc[2] - 2 * e2mc[3], e2mc[2] + 2 * e2mc[3]]]
        for ii in range(2):
            if tool_box.check_in(m_r[ii]):
                mc_title[ii] = ''
            else:
                mc_title[ii] = "_m" + str(ii+1)
            if tool_box.check_in(c_r[ii]):
                mc_title[ii + 2] = ''
            else:
                mc_title[ii + 2] = "_c" + str(ii+1)
        pic_mc = "".join(mc_title)

        pic_path = total_path + "result/cuts/sym/" + cut + "/" + str(cut_s) + pic_mc + ".eps"
        tool_box.mcplot(fg1, arr[0:3,:], fg2, arr[3:6,:], e1mc, e2mc, str(cut_s), 'max', pic_path)
        pic_path = total_path + "result/cuts/sym/" + cut + "/" + str(cut_s) + pic_mc + ".png"
        tool_box.mcplot(fg1, arr[0:3, :], fg2, arr[3:6, :], e1mc, e2mc, str(cut_s), 'max', pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = total_path + "result/cuts/sym/" + cut + "/total.npz"
    numpy.savez(mc_path, mc1, mc2)
    # mc1 = numpy.load(mc_path)['arr_0']
    # mc2 = numpy.load(mc_path)['arr_1']

    x1 = select[cut][1][0]
    x2 = select[cut][1][-1]

    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(121)
    ax1.errorbar(select[cut][1], mc1[0] - 1, mc1[1], c='coral', capsize=3, label='m1')
    ax1.errorbar(select[cut][1], mc2[0] - 1, mc2[1], c='royalblue', capsize=3, label='m2')
    ax1.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='lawngreen')
    ax1.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax1.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax1.set_xlabel("Cutoff")
    ax1.legend()
    # ax1.set_xscale('log')
    ax1.set_ylabel("m")

    ax2 = fig.add_subplot(122)
    ax2.errorbar(select[cut][1], mc1[2], mc1[3], c='coral', capsize=3, label='c1')
    ax2.errorbar(select[cut][1], mc2[2], mc2[3], c='royalblue', capsize=3, label='c2')
    ax2.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='lawngreen')
    ax2.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax2.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax2.set_xlabel("Cutoff")
    ax2.legend()
    # ax2.set_xscale('log')
    ax2.set_ylabel("c")

    namep = total_path + "result/cuts/sym/" + cut + "/total.eps"
    plt.savefig(namep)
    plt.close()

t2 = time.clock()
if rank == 0:
    print(t2-t1)

