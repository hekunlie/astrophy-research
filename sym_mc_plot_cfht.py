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
import shelve


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cut = argv[1]

g1num = 11
g2num = cpus
g1 = numpy.linspace(-0.0045, 0.0045, g1num)
g2 = numpy.linspace(-0.0095, 0.0095, g2num)
dg1 = g1[1]-g1[0]
dg2 = g2[1]-g2[0]

t1 = time.clock()
with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data" in path:
        total_path = path.split("=")[1]
    elif "cfht_res" in path:
        result_path = path.split("=")[1]

data_cache = result_path + "data_cache.npz"

data = numpy.load(data_cache)['arr_0']

fq = Fourier_Quad(48, 123)

# fg1 = data[:, 14]
# t_1 = numpy.isnan(fg1)
# idx_t1 = t_1 == False
#
# fg2 = data[:, 15]
# t_2 = numpy.isnan(fg2)
# idx_t2 = t_2 == False

peak = data[:, 4]
flux = data[:, 5]
hflux = data[:, 6]
area = data[:, 7]
harea = data[:, 8]
fsnr = data[:, 10]
fsnr_f = data[:, 11]
fg1 = data[:, 14]
fg2 = data[:, 15]

mg1 = data[:, 16]
mg2 = data[:, 17]
mn = data[:, 18]
mu = data[:, 19]
mv = data[:, 20]

cuts_num = 20

d_sort = numpy.sort(flux)
step = int(len(d_sort)/cuts_num)
fcut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(peak)
step = int(len(d_sort)/cuts_num)
pcut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(fsnr)
step = int(len(d_sort)/cuts_num)
fsnrcut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(fsnr_f)
step = int(len(d_sort)/cuts_num)
fsnrfcut = [d_sort[i*step] for i in range(cuts_num)]

select = {"flux":   (flux, fcut),     "peak":      (peak, pcut),
          "fsnr":   (fsnr, fsnrcut),  "fsnr_f":   (fsnr_f, fsnrfcut)}

if rank < g1num:
    idxg11 = fg1 >= g1[rank] - dg1/2
    idxg12 = fg1 <= g1[rank] + dg1/2

idxg21 = fg2 >= g2[rank] - dg2/2
idxg22 = fg2 <= g2[rank] + dg2/2

res_list = []
for tag, cut_s in enumerate(select[cut][1]):
    idx = select[cut][0] >= cut_s
    if rank < g1num:
        num1 = len(mg1[idx&idxg11&idxg12])
        g1_h, g1_sig = fq.fmin_g(mg1[idx&idxg11&idxg12], mn[idx&idxg11&idxg12], mu[idx&idxg11&idxg12], mode=1, bin_num=8)
    else:
        g1_h, g1_sig, num1 = -1, -1, -1

    num2 = len(mg2[idx&idxg21&idxg22])
    g2_h, g2_sig = fq.fmin_g(mg2[idx&idxg21&idxg22], mn[idx&idxg21&idxg22], mu[idx&idxg21&idxg22], mode=2, bin_num=8)
    res_list.append([g1_h, g1_sig, num1, g2_h, g2_sig, num2])

coll_res = comm.gather(res_list, root=0)

if rank == 0:
    # cache = shelve.open("/mnt/ddnfs/data_users/hkli/CFHT/result/cuts/peak/cache")
    # # cache['cache'] = coll_res
    # coll_res = cache['cache']
    # cache.close()
    mc1 = []
    mc2 = []
    for tag, cut_s in enumerate(select[cut][1]):
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
        xylim = (-0.008, 0.008, -0.012, 0.012)
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
    ax1.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='green')
    ax1.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax1.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax1.set_xlabel("Cutoff")
    ax1.legend()
    # ax1.set_xscale('log')
    ax1.set_ylabel("m")

    ax2 = fig.add_subplot(122)
    ax2.errorbar(x_coord, mc1[2], mc1[3], c='coral', capsize=3, label='c1')
    ax2.errorbar(x_coord, mc2[2], mc2[3], c='royalblue', capsize=3, label='c2')
    ax2.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='green')
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