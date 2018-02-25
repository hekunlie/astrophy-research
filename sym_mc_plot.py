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
import logging


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cut = argv[1]

logger = logging.getLogger()
logger.setLevel(logging.INFO)
logfile = '/lmc/selection_bias/logs/sym_%d_log.dat' %rank

lf = logging.FileHandler(logfile, 'w')
form = logging.Formatter('%(asctime)s - %(message)s')
lf.setFormatter(form)
logger.addHandler(lf)

t1 = time.clock()

shear = numpy.load("/lmc/selection_bias/parameters/shear.npz")
fg1 = shear["arr_0"]
fg2 = shear["arr_1"]

h5path = "/lmc/selection_bias/result/data/data_%d.hdf5"%rank
f = h5py.File(h5path, "r")
data = f["/data"].value
f.close()
logger.info('open data file data_%d.hdf5'%rank)
fq = Fourier_Quad(84, 152356)
nsig = 380.64

osnr = data[:, 7]
osnrcut = [0, 5, 10, 15, 20, 25, 35, 55, 75, 95, 110]
#
flux = data[:, 8]/nsig
fcut = [0, 20, 60, 100, 150, 200, 250, 300, 350, 400, 450]

peak = data[:, 9]/nsig
pcut = [0, 3, 4, 4.5, 5.5, 6.5, 7.5, 8.5, 10, 12, 15]
#
fsnr = data[:, 10]
fsnrcut = [0, 2.5, 3.5, 4.8, 6, 7, 8, 10, 15, 20, 25]

fsnr4 = data[:, 11]
fsnr4cut = [0, 2.5, 3.5, 4.8, 6, 7, 8, 10, 15, 20, 25]

fsnr_c = data[:, 12]
fsnrccut = [0, 2.5, 3.5, 4.8, 6, 7, 8, 10, 15, 20, 25]
#
fsnr_c4 = data[:, 13]
fsnrc4cut = [0, 2.5, 3.5, 4.8, 6, 7, 8, 10, 15, 20, 25]

snr = data[:, 14]
snrcut = [0, 5, 10, 15, 20, 25, 35, 55, 75, 95, 110]

sesnr = data[:, 15]
secut = [0, 13, 15, 18, 21, 25, 30, 40, 60, 90, 120]

select = {"sesnr": (sesnr, secut), 'osnr':(osnr, osnrcut),"fsnr": (fsnr, fsnrcut), "flux": (flux, fcut), "peak": (peak, pcut), "fsnr4": (fsnr4, fsnr4cut),
          "fsnr_c": (fsnr_c, fsnrccut), "fsnr_c4": (fsnr_c4, fsnrc4cut), 'snr':(snr, snrcut)}

res_arr = numpy.zeros((6, len(select[cut][1])))

mg1 = data[:, 2]
mg2 = data[:, 3]
mn = data[:, 4]
mu = data[:, 5]
mv = data[:, 6]
logger.info('begin to do cutoff')
for tag, cut_s in enumerate(select[cut][1]):
    t_c1 = time.clock()
    idx = select[cut][0] > cut_s
    num = len(mg1[idx])
    logger.info("num: %d, %.2f"%(num, cut_s))
    g1_h, g1_sig = fq.fmin_g(mg1[idx], mn[idx], mu[idx], mode=1, bin_num=8)
    logger.info('g1')
    g2_h, g2_sig = fq.fmin_g(mg2[idx], mn[idx], mu[idx], mode=2, bin_num=8)
    logger.info('g2')
    res_arr[:, tag] = numpy.array([g1_h, g1_sig, num, g2_h, g2_sig, num])
    t_c2 = time.clock()
    logger.info("Procs: %d, %d, cuts: %.2f, time: %.2f"%(rank, tag, cut_s, t_c2-t_c1))

if rank > 0:
    # pass
    comm.Send(res_arr, dest=0, tag=rank)
    logger.info('%d sent information'%rank)
else:
    for procs in range(1, cpus):
        logger.info('begin receive from %d'%procs)
        recvs = numpy.empty((6, len(select[cut][1])), dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        res_arr = numpy.column_stack((res_arr, recvs))
        logger.info('received from %d' % procs)
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
        data_path = "/lmc/selection_bias/result/cuts/sym/" + cut + "/" + str(cut_s)+".npz"
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

        pic_path = "/lmc/selection_bias/result/cuts/sym/" + cut + "/" + str(cut_s) + pic_mc + ".eps"
        tool_box.mcplot(fg1, arr[0:3,:], fg2, arr[3:6,:], e1mc, e2mc, str(cut_s), 'max', pic_path)
        pic_path = "/lmc/selection_bias/result/cuts/sym/" + cut + "/" + str(cut_s) + pic_mc + ".png"
        tool_box.mcplot(fg1, arr[0:3, :], fg2, arr[3:6, :], e1mc, e2mc, str(cut_s), 'max', pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = "/lmc/selection_bias/result/cuts/sym/" + cut + "/total.npz"
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

    namep = "/lmc/selection_bias/result/cuts/sym/" + cut + "/total.eps"
    plt.savefig(namep)
    plt.close()

t2 = time.clock()
if rank == 0:
    print(t2-t1)

