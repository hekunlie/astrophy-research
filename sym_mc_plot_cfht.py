import matplotlib
matplotlib.use("Agg")
import numpy
import matplotlib.pyplot as plt
import os
from sys import path
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import time
from Fourier_Quad import Fourier_Quad
from sys import argv
from mpi4py import MPI
import tool_box
import shelve


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cut = argv[1]

t1 = time.clock()
with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        total_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_field_path" in path:
        field_path = path.split("=")[1]
    elif "cfht_cut_path" in path:
        cut_path = path.split("=")[1]

g_bin_path = result_path + "g_bin.npz"
g_data = numpy.load(g_bin_path)
g1 = g_data['arr_0']
g2 = g_data['arr_1']
g1num = len(g1)
g2num = len(g2)
# the length of the interval
dg1 = g1[1] - g1[0]
dg2 = g2[1] - g2[0]

g_trues = [g1, g2]
dgs = [dg1, dg2]
gnums = [g1num, g2num]

if cpus < max(gnums):
    if rank == 0:
        print("Number of threads (%d) is smaller than the field bin number (%d)!"%(cpus, max(gnums)))
    exit()

select_save = shelve.open(cut_path+"cut_dict")
select = select_save['dict']
select_save.close()
if cut not in select.keys():
    if rank == 0:
        print("%s is not in the cutoff dict!"%cut)
    exit()
else:
    cuts_num = len(select[cut])

fq = Fourier_Quad(48, 123)
sp = (cuts_num, 6)
cut_result = numpy.zeros(sp)
for i in range(2):
    g_true = g_trues[i]
    dg = dgs[i]
    data_cache = result_path + "data/g%d_%d.npz" % (i + 1, rank)
    if os.path.exists(data_cache):
        data = numpy.load(data_cache)['arr_0']

        peak = data[:, 4]  # &bi_idx&field_idx]
        flux = data[:, 5]
        hflux = data[:, 6]
        area = data[:, 7]
        harea = data[:, 8]
        flux2 = data[:, 10]
        flux_alt = data[:, 11]
        field_g1 = data[:, 14]
        field_g2 = data[:, 15]
        MG1 = data[:, 16]
        MG2 = data[:, 17]
        MN = data[:, 18]
        MU = data[:, 19]
        # MV = data[:, 20][idx & idxa]
        DE1 = MN + MU
        DE2 = MN - MU

        MGs = [MG1, MG2]
        DEs = [DE1, DE2]
        select_cri = {"flux": flux, "hflux": hflux, "area": area, "harea": harea, "peak": peak,
                      "flux2": flux2, "flux_alt": flux_alt}
        for tag, cut_s in enumerate(select[cut]):
            idx_c = select_cri[cut] >= cut_s
            mg = MGs[i][idx_c]
            de = DEs[i][idx_c]
            gal_num = len(mg)
            g_h, g_sig = fq.fmin_g(mg, de, bin_num=12)
            cut_result[tag, i*3:i*3+3] = g_h, g_sig, gal_num

if rank > 0:
    comm.Send(cut_result, dest=0, tag=rank)
else:
    res_list = [cut_result]
    for procs in range(1, cpus):
        recvs = numpy.empty(sp, dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        res_list.append(recvs)
    # cache = shelve.open("/mnt/ddnfs/data_users/hkli/CFHT/result/cuts/peak/cache")
    # # cache['cache'] = coll_res
    # coll_res = cache['cache']
    # cache.close()
    mc1 = []
    mc2 = []
    mcs = [[], []]
    for tag, cut_s in enumerate(select[cut]):
        arr = []
        for i in range(max(gnums)):
            arr.append(res_list[i][tag])
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
        data_path = result_path + "cuts/" + cut + "/" + str(round(cut_s, 4)) + ".npz"
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
        xylim = (-0.008, 0.008, -0.009, 0.009)
        pic_path = result_path + "cuts/" + cut + "/" + str(round(cut_s, 4)) + pic_mc + ".eps"
        tool_box.mcplot(g1, arr[0:3, 0:g1num], g2, arr[3:6, 0:g2num], e1mc, e2mc, str(round(cut_s, 4)), 'max', xylim,pic_path)
        pic_path = result_path + "cuts/" + cut + "/" + str(round(cut_s, 4)) + pic_mc + ".png"
        tool_box.mcplot(g1, arr[0:3, 0:g1num], g2, arr[3:6, 0:g2num], e1mc, e2mc, str(round(cut_s, 4)), 'max', xylim,pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = result_path + "cuts/" + cut + "/total.npz"
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

    namep = result_path + "cuts/"+ cut + "/total.eps"
    plt.savefig(namep)
    plt.close()

t2 = time.clock()
if rank == 0:
    print(t2 - t1)
