# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import os
from sys import path
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
from Fourier_Quad import Fourier_Quad
from sys import argv
import numpy
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

del_bin, bin_num, cho, cho_thre = argv[1],argv[2],argv[3],argv[4]

ts = time.clock()

del_bin = int(del_bin)
cho_thre = float(cho_thre)
bin_num = int(bin_num)

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_pic_path" in path:
        pic_path = path.split("=")[1]

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

measure_res = []
for i in range(2):
    g_true = g_trues[i]
    dg = dgs[i]
    data_cache = result_path + "data/g%d_%d.npz"%(i+1,rank)
    if os.path.exists(data_cache):
        data = numpy.load(data_cache)['arr_0']
        stars = data[:, 3]
        idx = stars >= 12

        peak = data[:, 4][idx]
        flux = data[:, 5][idx]
        hflux = data[:, 6][idx]
        area = data[:, 7][idx]
        harea = data[:, 8][idx]
        flux2 = data[:, 10][idx]
        flux_alt = data[:, 11][idx]
        field_g1 = data[:, 14][idx]
        field_g2 = data[:, 15][idx]
        MG1 = data[:, 16][idx]
        MG2 = data[:, 17][idx]
        MN = data[:, 18][idx]
        MU = data[:, 19][idx]
        MV = data[:, 20][idx]
        # be careful that the "MU" defined in FRESH is the different from that in ours
        # MN + MU for our definition of MU and MV which is the same as those in the paper Zhang et al. 2017 ApJ, 834:8
        DE1 = MN - MU
        DE2 = MN + MU

        selects = {"peak": peak, "flux2": flux2, "flux_alt": flux_alt, "flux": flux}
        sel_idx = selects[cho] >= cho_thre

        MGs = [MG1, MG2]
        DEs = [DE1, DE2]
        fgs = [field_g1, field_g2]

        # fg = fgs[i]
        # mg = MGs[i]
        # de = DEs[i]

        # idx1 = fgs[i] >= g_true[rank] - dg/2
        # idx2 = fgs[i] <= g_true[rank] + dg/2

        mg = MGs[i][sel_idx]
        de = DEs[i][sel_idx]

        pic = pic_path + "%s_%d_%.2f_g%d_%d.png"%(cho,del_bin, cho_thre, i+1,rank)
        estg, sig = Fourier_Quad(48,123).fmin_g_new(g=mg, nu=de, bin_num=bin_num, ig_num=del_bin, pic_path=pic)
        gal_num = len(mg)
        field_g = g_true[rank]

    else:
        estg, sig, gal_num, field_g = -1,-1,-1, -1
    measure_res.extend([estg, sig, gal_num, field_g])

gather_data = comm.gather(measure_res, root=0)

if rank == 0:
    res_data = numpy.array(gather_data)
    mcs = []
    for i in range(2):
        mg = res_data[0:gnums[i], i*3+0+i]
        sig = res_data[0:gnums[i], i*3+1+i]
        num = res_data[0:gnums[i], i*3+2+i]
        mc = tool_box.data_fit(g_trues[i], mg, sig)
        mcs.append(mc)
        for p in range(gnums[i]):
            print("num: %4.1f W, g%d: %8.5f, m_g: %8.5f, sig: %8.5f, devi: %4.2f * e^-4, shape noise: %6.4f"
                  %(num[p]/10000, i+1, g_trues[i][p], mg[p], sig[p], 10000*(mg[p]-g_trues[i][p]), numpy.sqrt(num[p])*sig[p]))
        print("\n")
    final_cache = result_path + "%d_%s_%.2f_final_cache.npz"%(del_bin, cho, cho_thre)
    final_cache = tool_box.file_name(final_cache)
    numpy.savez(final_cache, numpy.array(mcs), res_data)
    e1mc = mcs[0]
    e2mc = mcs[1]
    if e1mc[0] - 1 - 2 * e1mc[1] < 0 < e1mc[0] - 1 + 2 * e1mc[1]:
        m1_b = "no m1 bias"
    else:
        m1_b = "m1 bias"
    if e2mc[0] - 1 - 2 * e2mc[1] < 0 < e2mc[0] - 1 + 2 * e2mc[1]:
        m2_b = "no m2 bias"
    else:
        m2_b = "m2 bias"

    if e1mc[2] - 2 * e1mc[3] < 0 < e1mc[2] + 2 * e1mc[3]:
        c1_b = "no c1 bias"
    else:
        c1_b = "c1 bias"
    if e2mc[2] - 2 * e2mc[3] < 0 < e2mc[2] + 2 * e2mc[3]:
        c2_b = "no c2 bias"
    else:
        c2_b = "c2 bias"
    print("%10s: %8.5f (%6.5f), %10s: %10.6f (%.6f)" % (m1_b, e1mc[0] - 1, e1mc[1], c1_b, e1mc[2], e1mc[3]))
    print("%10s: %8.5f (%6.5f), %10s: %10.6f (%.6f)" % (m2_b, e2mc[0] - 1, e2mc[1], c2_b, e2mc[2], e2mc[3]))

    nm = pic_path + "%d_%s_%.2f_mc.png"%(del_bin,cho, cho_thre)
    nm = tool_box.file_name(nm)
    g1_bound = 1.4*numpy.abs(g1).max()
    g2_bound = 1.4*numpy.abs(g2).max()
    xy_lim = [-g1_bound, g1_bound, -g2_bound, g2_bound]
    tool_box.mcplot(g1, res_data.T[0:3,0:g1num], g2, res_data.T[4:7, 0:g2num], e1mc, e2mc, str(cho_thre), 'max', xy_lim, nm)
te = time.clock()
if rank == 0:
    print(te-ts)


