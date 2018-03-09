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
for i in range(5):
    h5path = total_path + "result/data/data_%d_%d.hdf5"%(rank,i)
    f = h5py.File(h5path, "r")
    if i == 0:
        data = f["/data"].value
    else:
        data = numpy.row_stack((data, f["/data"].value))
    f.close()

fq = Fourier_Quad(90, 152356)
nsig = 380.64
method = 'mean'
cuts_num = 20
osnr = data[:, 7]
detect_idx = osnr >0
d_sort = numpy.sort(osnr[detect_idx])
step = int(len(d_sort)/cuts_num)
osnrcut = [d_sort[i*step] for i in range(cuts_num)]

flux = data[:, 8]/nsig
d_sort = numpy.sort(flux[detect_idx])
step = int(len(d_sort)/cuts_num)
fcut = [d_sort[i*step] for i in range(cuts_num)]

peak = data[:, 9]/nsig
d_sort = numpy.sort(peak[detect_idx])
step = int(len(d_sort)/cuts_num)
pcut = [d_sort[i*step] for i in range(cuts_num)]

fsnr_c = data[:, 10]
d_sort = numpy.sort(fsnr_c[detect_idx])
step = int(len(d_sort)/cuts_num)
fsnrccut = [d_sort[i*step] for i in range(cuts_num)]

snr = data[:, 11]
d_sort = numpy.sort(snr[detect_idx])
step = int(len(d_sort)/cuts_num)
snrcut = [d_sort[i*step] for i in range(cuts_num)]

fpath = total_path + "result/data/fsnr_%d.npz"%rank
f = numpy.load(fpath)["arr_0"]
fsnr = f[:, 0]
d_sort = numpy.sort(fsnr[detect_idx])
step = int(len(d_sort)/cuts_num)
fsnrcut = [d_sort[i*step] for i in range(cuts_num)]

fsnr_f = f[:, 1]
d_sort = numpy.sort(fsnr_f[detect_idx])
step = int(len(d_sort)/cuts_num)
fsnrfcut = [d_sort[i*step] for i in range(cuts_num)]

spath = total_path + "result/data/sex_%d.npz"%rank
sesnr = numpy.load(spath)["arr_0"][:,0]
d_sort = numpy.sort(sesnr[sesnr>0])
step = int(len(d_sort)/cuts_num)
secut = [d_sort[i*step] for i in range(cuts_num)]


select = {"sesnr": (sesnr, secut), 'osnr':(osnr, osnrcut),
          "fsnr": (fsnr, fsnrcut), "flux": (flux, fcut),
          "peak": (peak, pcut), "fsnr_c": (fsnr_c, fsnrccut),
          "fsnr_f": (fsnr_f, fsnrfcut), 'snr':(snr, snrcut)}

res_arr = numpy.zeros((6, len(select[cut][1])))

mg1 = data[:, 2]
mg2 = data[:, 3]
mn = data[:, 4]
mu = data[:, 5]
mv = data[:, 6]

for tag, cut_s in enumerate(select[cut][1]):
    idx = select[cut][0] > cut_s
    num = len(mg1[idx])

    # g1_h, g1_sig = fq.fmin_g(mg1[idx], mn[idx], mu[idx], mode=1, bin_num=8)
    # g2_h, g2_sig = fq.fmin_g(mg2[idx], mn[idx], mu[idx], mode=2, bin_num=8)

    weight = 1#select[cut][0][idx]
    g1_h = numpy.mean(mg1[idx] / weight) / numpy.mean(mn[idx] / weight)
    g1_sig = numpy.sqrt(numpy.mean((mg1[idx] / weight) ** 2) / (numpy.mean(mn[idx] / weight)) ** 2) / numpy.sqrt(num)

    g2_h = numpy.mean(mg2[idx] / weight) / numpy.mean(mn[idx] / weight)
    g2_sig = numpy.sqrt(numpy.mean((mg2[idx] / weight) ** 2) / (numpy.mean(mn[idx] / weight)) ** 2) / numpy.sqrt(num)

    res_arr[:, tag] = numpy.array([g1_h, g1_sig, num, g2_h, g2_sig, num])
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
        data_path = total_path + "result/cuts/%s/"%method + cut + "/" + str(round(cut_s,4))+".npz"
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

        pic_path = total_path + "result/cuts/%s/"%method + cut + "/" + str(round(cut_s,4)) + pic_mc + ".eps"
        tool_box.mcplot(fg1, arr[0:3,:], fg2, arr[3:6,:], e1mc, e2mc, str(round(cut_s,4)), 'max', pic_path)
        pic_path = total_path + "result/cuts/%s/"%method + cut + "/" + str(round(cut_s,4)) + pic_mc + ".png"
        tool_box.mcplot(fg1, arr[0:3, :], fg2, arr[3:6, :], e1mc, e2mc, str(round(cut_s,4)), 'max', pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = total_path + "result/cuts/%s/"%method + cut + "/total.npz"
    numpy.savez(mc_path, mc1, mc2)
    # mc1 = numpy.load(mc_path)['arr_0']
    # mc2 = numpy.load(mc_path)['arr_1']

    x1 = 0
    x2 = 1
    x_coord = [i * 1/cuts_num for i in range(cuts_num)]

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

    namep = total_path + "result/cuts/%s/"%method + cut + "/total.eps"
    plt.savefig(namep)
    plt.close()

t2 = time.clock()
if rank == 0:
    print(t2-t1)

