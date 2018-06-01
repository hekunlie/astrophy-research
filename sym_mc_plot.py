import matplotlib
matplotlib.use("Agg")
import numpy
import matplotlib.pyplot as plt
import os
from sys import path,argv
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import time
from Fourier_Quad import Fourier_Quad
from mpi4py import MPI
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cut = argv[1]

g1num = cpus-6
g2num = cpus
g1 = numpy.linspace(-0.004, 0.004, g1num)
g2 = numpy.linspace(-0.0055, 0.0055, g2num)
dg1 = g1[1]-g1[0]
dg2 = g2[1]-g2[0]

t1 = time.clock()
with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_pic_path" in path:
        pic_path = path.split("=")[1]
    elif "cfht_cut_path" in path:
        cut_path = path.split("=")[1]

fq = Fourier_Quad(48, 123)

g1_data_path = result_path + "g1_%d.npz"%rank
g2_data_path = result_path + "g2_%d.npz"%rank

g1num = cpus-6
g2num = cpus
g1 = numpy.linspace(-0.004, 0.004, g1num)
g2 = numpy.linspace(-0.0055, 0.0055, g2num)

cuts_num = 20
star_thre = 16

cut_off_scale = numpy.load(result_path+cut+'.npz')["arr_0"]

if rank < g1num:
    g1_data = numpy.load(g1_data_path)['arr_0']
    n_star_1 = g1_data[:, 3]
    idx_1 = n_star_1 >= star_thre
    nsig_1 = g1_data[:, 4][idx_1]
    flux_1 = g1_data[:, 5][idx_1]
    fsnr_1 = numpy.sqrt(g1_data[:, 10][idx_1])
    FG1_1 = g1_data[:, 16][idx_1]
    FN_1 = g1_data[:, 18][idx_1]
    FU_1 = g1_data[:, 19][idx_1]
    snr08_1 = numpy.sqrt(flux_1)/nsig_1

    select_1 = {"flux": flux_1, "snr": snr08_1, "fsnr": fsnr_1}

g2_data = numpy.load(g2_data_path)['arr_0']
n_star_2 = g2_data[:, 3]
idx_2 = n_star_2 >= star_thre

nsig_2 = g2_data[:, 4][idx_2]
flux_2 = g2_data[:, 5][idx_2]
fsnr_2 = numpy.sqrt(g2_data[:, 10][idx_2])
FG2_2 = g2_data[:, 17][idx_2]
FN_2 = g2_data[:, 18][idx_2]
FU_2 = g2_data[:, 19][idx_2]
snr08_2 = numpy.sqrt(flux_2)/nsig_2
select_2 = {"flux": flux_2, "snr": snr08_2, "fsnr": fsnr_2}

res_list = []
for i in range(cuts_num):
    if rank < g1num:
        idx = select_1[cut] >= cut_off_scale[i]
        num1 = len(FG1_1[idx])
        g1_h, g1_sig = fq.fmin_g(FG1_1[idx], FN_1[idx], FU_1[idx], mode=1, bin_num=8)
    else:
        g1_h, g1_sig, num1 = -1, -1, -1

    idx = select_2[cut] >= cut_off_scale[i]
    num2 = len(FG2_2[idx])
    g2_h, g2_sig = fq.fmin_g(FG2_2[idx], FN_2[idx], FU_2[idx], mode=2, bin_num=8)
    res_list.append([g1_h, g1_sig, num1, g2_h, g2_sig, num2])

coll_res = comm.gather(res_list, root=0)

if rank == 0:
    mc1 = []
    mc2 = []
    for tag, cut_s in enumerate(cut_off_scale):
        arr = []
        for i in range(cpus):
            arr.append(coll_res[i][tag])
        arr = numpy.array(arr).T
        y1_data = arr[0, 0:g1num]
        y1_err = arr[1, 0:g1num]
        y2_data = arr[3, 0:g2num]
        y2_err = arr[4, 0:g2num]

        e1mc = tool_box.data_fit(g1, y1_data, y1_err)
        mc1.append(e1mc)
        e2mc = tool_box.data_fit(g2, y2_data, y2_err)
        mc2.append(e2mc)

        mc = numpy.array([e1mc, e2mc])
        data_path = cut_path + "%s/%.4f.npz"%(cut,cut_s)

        numpy.savez(data_path, arr, mc)

        mc_title = ['0', '0', '0', '0']
        m_r = [[e1mc[0] - 1 - 2 * e1mc[1], e1mc[0] - 1 + 2 * e1mc[1]],
               [e2mc[0] - 1 - 2 * e2mc[1], e2mc[0] - 1 + 2 * e2mc[1]]]
        c_r = [[e1mc[2] - 2 * e1mc[3], e1mc[2] + 2 * e1mc[3]], [e2mc[2] - 2 * e2mc[3], e2mc[2] + 2 * e2mc[3]]]
        for ii in range(2):
            if tool_box.check_in(m_r[ii]):
                mc_title[ii] = ''
            else:
                mc_title[ii] = "_m" + str(ii + 1)
            if tool_box.check_in(c_r[ii]):
                mc_title[ii + 2] = ''
            else:
                mc_title[ii + 2] = "_c" + str(ii + 1)
        pic_mc = "".join(mc_title)
        xylim = (-0.0065, 0.0065, -0.0075, 0.0075)
        pic_path = cut_path + cut + "/" + str(round(cut_s, 4)) + pic_mc + ".eps"
        tool_box.mcplot(g1, arr[0:3, 0:g1num], g2, arr[3:6, 0:g2num], e1mc, e2mc, str(round(cut_s, 4)), 'max', xylim,pic_path)
        pic_path = cut_path + cut + "/" + str(round(cut_s, 4)) + pic_mc + ".png"
        tool_box.mcplot(g1, arr[0:3, 0:g1num], g2, arr[3:6, 0:g2num], e1mc, e2mc, str(round(cut_s, 4)), 'max', xylim,pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = cut_path + cut + "/total.npz"
    numpy.savez(mc_path, mc1, mc2)
    # mc1 = numpy.load(mc_path)['arr_0']
    # mc2 = numpy.load(mc_path)['arr_1']

    x1 = 0
    x2 = 1
    x_coord = [i * 1. / cuts_num for i in range(cuts_num)]

    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(121)
    ax1.errorbar(x_coord, mc1[0] - 1, mc1[1], c='coral', capsize=3, label='m1')
    ax1.errorbar(x_coord, mc2[0] - 1, mc2[1], c='royalblue', capsize=3, label='m2')
    ax1.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='grey')
    ax1.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax1.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax1.set_xlabel("Cutoff")
    ax1.legend()
    # ax1.set_xscale('log')
    ax1.set_ylabel("m")

    ax2 = fig.add_subplot(122)
    ax2.errorbar(x_coord, mc1[2], mc1[3], c='coral', capsize=3, label='c1')
    ax2.errorbar(x_coord, mc2[2], mc2[3], c='royalblue', capsize=3, label='c2')
    ax2.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='grey')
    ax2.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax2.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax2.set_xlabel("Cutoff")
    ax2.legend()
    # ax2.set_xscale('log')
    ax2.set_ylabel("c")

    namep = cut_path + cut + "/total.eps"
    plt.savefig(namep)
    plt.close()

t2 = time.clock()
if rank == 0:
    print(t2 - t1)
