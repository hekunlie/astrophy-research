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

data_cache = result_path + "data_cache.npz"
data = numpy.load(data_cache)['arr_0']

if rank == 0:
    print("Totally, %d galaxies are detected"%len(data))


fg1_max, fg1_min = numpy.max(data[:, 14]),numpy.min(data[:, 14])
fg2_max, fg2_min = numpy.max(data[:, 15]),numpy.min(data[:, 15])

# binary_tag = binary_data[:, 0]
# field_lab = binary_data[:, 1]
# expo_lab = binary_data[:, 2]
# chip_lab = binary_data[:, 3]
# '1' means binary or triple
# bi_idx = binary_tag != 1
# exclude some fields
# field_idx = field_lab != 41100
# if rank == 0:
#     print("Binary_detect", len(binary_tag) - len(binary_tag[bi_idx]))
#     print("Field excluded contains:", len(field_lab) - len(field_lab[field_idx]))

n_star = data[:, 3]
idx = n_star >= 20
area = data[:, 7]
idxa = area >= 5

peak = data[:, 4][idx&idxa]#&bi_idx&field_idx]
flux = data[:, 5][idx&idxa]
hflux = data[:, 6][idx&idxa]
area = data[:, 7][idx&idxa]
harea = data[:, 8][idx&idxa]
fsnr = data[:, 10][idx&idxa]
fsnr_f = data[:, 11][idx&idxa]
fg1 = data[:, 14][idx&idxa]
fg2 = data[:, 15][idx&idxa]
FG1 = data[:, 16][idx&idxa]
FG2 = data[:, 17][idx&idxa]
FN = data[:, 18][idx&idxa]
FU = data[:, 19][idx&idxa]
FV = data[:, 20][idx&idxa]
DE1 = FN - FU
DE2 = FN + FU
selects = {"peak": peak, "fsnr": fsnr, "fsnr_f": fsnr_f, "flux": flux}
sel_idx = selects[cho] >= cho_thre

g1num = cpus
g2num = cpus
g1 = numpy.linspace(-0.005, 0.005, g1num)
g2 = numpy.linspace(-0.005, 0.005, g2num)

# the length of the interval
dg1 = g1[1] - g1[0]
dg2 = g2[1] - g2[0]

if rank == 0:
    print("nstar: %d ~ %d\n"%(numpy.min(n_star), numpy.max(n_star)))
    print("fg1: %.4f ~ %.4f (%.4f ~ %.4f)\n"%(fg1_min, fg1_max, numpy.min(fg1), numpy.max(fg1)))
    print("fg2: %.4f ~ %.4f (%.4f ~ %.4f)\n" % (fg2_min, fg2_max, numpy.min(fg2), numpy.max(fg2)))
    print("peak: %.2f ~ %.2f\n" % (numpy.min(peak), numpy.max(peak)))
    print("flux: %.2f ~ %.2f\n" % (numpy.min(flux), numpy.max(flux)))
    print("area: %.2f ~ %.2f\n" % (numpy.min(area), numpy.max(area)))

if rank < g1num:
    idx11 = fg1 >= g1[rank] - dg1/2
    idx12 = fg1 <= g1[rank] + dg1/2

    mg1 = FG1[idx11&idx12&sel_idx]
    de1 = DE1[idx11&idx12&sel_idx]

    pic_1 = pic_path + "%d_%s_%.2f_g1_%d.png"%(del_bin, cho, cho_thre, rank)
    esg1, sig1 = Fourier_Quad(48,123).fmin_g_new(g=mg1, nu=de1, mode=1, bin_num=bin_num, ig_num=del_bin, pic_path=pic_1)
    field_g1 = g1[rank]
    num1 = len(mg1)
else:
    esg1, sig1, num1, field_g1 = -1,-1,-1,-1

idx21 = fg2 >= g2[rank] - dg2/2
idx22 = fg2 <= g2[rank] + dg2/2

mg2 = FG2[idx21&idx22&sel_idx]
de2 = DE2[idx21&idx22&sel_idx]

pic_2 = pic_path + "%d_%s_%.2f_g2_%d.png"%(del_bin, cho, cho_thre, rank)
esg2, sig2 = Fourier_Quad(48,123).fmin_g_new(g=mg2,nu=de2, mode=2,bin_num=bin_num,ig_num=del_bin,pic_path=pic_2)
field_g2 = g2[rank]
num2 = len(mg2)

send_data = [esg1, sig1, num1, esg2, sig2, num2, field_g1, field_g2]
gather_data = comm.gather(send_data, root=0)

if rank == 0:
    g_data = numpy.array(gather_data)
    g = [g1,g2]
    mcs = []
    gnum = [g1num, g2num]
    for i in range(2):
        mg = g_data[0:gnum[i], i*3+0]
        sig = g_data[0:gnum[i], i*3+1]
        num = g_data[0:gnum[i], i*3+2]
        mc = tool_box.data_fit(g[i], mg, sig)
        mcs.append(mc)
        for p in range(gnum[i]):
            print("num: %4.1f W, g1: %8.5f, m_g: %8.5f, sig: %8.5f, devi: %4.2f * e^-4, shape noise: %6.4f"
                  %(num[p]/10000, g[i][p], mg[p], sig[p], 10000*(mg[p]-g[i][p]), numpy.sqrt(num[p])*sig[p]))
        print("\n")
    final_cache = result_path + "%d_%s_%.2f_final_cache.npz"%(del_bin, cho, cho_thre)
    final_cache = tool_box.file_name(final_cache)
    numpy.savez(final_cache, numpy.array(mcs), g_data)
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
    xy_lim = [-0.008, 0.008, -0.010, 0.010]
    tool_box.mcplot(g1, g_data.T[0:3,0:g1num], g2, g_data.T[3:6, 0:g2num], e1mc, e2mc, str(cho_thre), 'max', xy_lim, nm)
te = time.clock()
if rank == 0:
    print(te-ts)


