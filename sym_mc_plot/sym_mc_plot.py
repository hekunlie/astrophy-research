import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import time
from Fourier_Quad import Fourier_Quad
from sys import argv
from mpi4py import MPI
import tool_box
import h5py
import numpy
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cut, filter_name, source, sig, r_thresh = argv[1], argv[2], argv[3], argv[4], float(argv[5])

t1 = time.clock()

env_path = "%s/work/envs/envs.dat"%my_home

path_items = tool_box.config(env_path, ['get','get','get'], [['selection_bias', "%s_path"%source, '1'],
                                                             ['selection_bias', "%s_path_result"%source, '1'],
                                                             ['selection_bias', "%s_path_para"%source, '1']])

total_path, result_path, para_path = path_items
shear_path = para_path + "shear.npz"
shear = numpy.load(shear_path)
fg1 = shear["arr_0"]
fg2 = shear["arr_1"]


sex_path = total_path + "result/data/%s/sex_%d.npz"%(filter_name,rank)
para_h5_path = para_path + "para_%d.hdf5"%rank
para_h5 = h5py.File(para_h5_path,"r")
mag_true = para_h5["/mag"].value
para_h5.close()
# rfactor_path = result_path + "data/resolution_factor_%d.npz"%rank

if rank == 0:
    log = "START -- %s, %7s, %9s,      , %5.2f, %12s"%(total_path.split("/")[-2], filter_name, cut, r_thresh, argv[0])
    logger = tool_box.get_logger("%s/work/selection_bias/sym_mc_plot/cutoff.dat"%my_home)
    logger.info(log)

shear_esti_h5path = result_path + "data/data_%d.hdf5"%rank
shear_esti_file = h5py.File(shear_esti_h5path, "r")
shear_esti_data = shear_esti_file["/data"].value
shear_esti_file.close()

fq = Fourier_Quad(60, 152356)
noise_sig = 60

MG1 = shear_esti_data[:, 2]
MG2 = shear_esti_data[:, 3]
MN = shear_esti_data[:, 4]
MU = shear_esti_data[:, 5]
MV = shear_esti_data[:, 6]
# be careful that the "MU" defined in FRESH is different from that in ours
# MN + (-) MU for our definition (g1(2)) of MU and MV which is the same as
# those in the paper Zhang et al. 2017 ApJ, 834:8
DE1 = MN + MU
DE2 = MN - MU

cuts_num = 10

fourier_idx = ["flux2", "flux_alt", "flux", "snr",
               "flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5",
               "flux_ex1", "flux_ex2", "flux_ex3", "flux_ex4", "flux_ex5",
               "snr_ex1", "snr_ex2", "snr_ex3", "snr_ex4", "snr_ex5"]
sex_idx = ["sex_snr", "sex_area", "snr_auto",  "mag_auto"]
if cut in fourier_idx:

    fq_snr_h5path = result_path + "data/data_%s/data_%d.hdf5" % (sig, rank)
    fq_snr_file = h5py.File(fq_snr_h5path, "r")
    fq_snr_data = fq_snr_file["/data"].value
    fq_snr_file.close()
    # it should be
    detected = fq_snr_data[:, 20] > -1
    if "snr_ex" in cut:
        flux_data = fq_snr_data[:, fourier_idx.index(cut)-5]
        area_data = fq_snr_data[:, fourier_idx.index(cut)]
        idx = area_data == 0
        area_data[idx] = 1
        cut_data = flux_data/numpy.sqrt(area_data)/noise_sig
    else:
        cut_data = fq_snr_data[:, fourier_idx.index(cut)]


elif cut in sex_idx:
    sex_data = numpy.load(sex_path)["arr_0"]
    detected = sex_data[:, 0] > 0
    if "snr_auto" in cut:
        flux_auto = sex_data[:, 1]
        flux_err = sex_data[:, 2]
        err_idx = flux_err == 0
        flux_err[err_idx] = 1
        cut_data = flux_auto / flux_err
    else:
        cut_data = sex_data[:, sex_idx.index(cut)]

cut_data_sort = numpy.sort(cut_data[detected])
cut_data_step = int(len(cut_data_sort) / cuts_num)
cutoff_scalar = [cut_data_sort[i * cut_data_step] for i in range(cuts_num)]
# gathering for check
cutoff_scalar_total = comm.gather(cutoff_scalar, root=0)
res_arr = numpy.zeros((6, len(cutoff_scalar)))

for tag, cut_s in enumerate(cutoff_scalar):
    idx = cut_data >= cut_s

    nm1 = MG1[detected&idx]
    de1 = DE1[detected&idx]
    nm2 = MG2[detected&idx]
    de2 = DE2[detected&idx]
    num = len(nm1)
    g1_h, g1_sig = fq.fmin_g_new(nm1, de1, bin_num=8)
    g2_h, g2_sig = fq.fmin_g_new(nm2, de2, bin_num=8)

    # weight = select[cut][0][idx]
    # g1_h = numpy.mean(mg1[idx] / weight) / numpy.mean(mn[idx] / weight)
    # g1_sig = numpy.sqrt(numpy.mean((mg1[idx] / weight) ** 2) / (numpy.mean(mn[idx] / weight)) ** 2) / numpy.sqrt(num)
    # g2_h = numpy.mean(mg2[idx] / weight) / numpy.mean(mn[idx] / weight)
    # g2_sig = numpy.sqrt(numpy.mean((mg2[idx] / weight) ** 2) / (numpy.mean(mn[idx] / weight)) ** 2) / numpy.sqrt(num)

    res_arr[:, tag] = numpy.array([g1_h, g1_sig, num, g2_h, g2_sig, num])

if rank > 0:
    comm.Send(res_arr, dest=0, tag=rank)
else:
    for procs in range(1, cpus):
        recvs = numpy.empty((6, len(cutoff_scalar)), dtype=numpy.float64)
        comm.Recv(recvs, source=procs, tag=procs)
        res_arr = numpy.column_stack((res_arr, recvs))

    numpy.savez(result_path+"cuts/cache_%s.npz"%cut, res_arr, cutoff_scalar_total)

    mc1 = []
    mc2 = []
    for tag, cut_s in enumerate(cutoff_scalar):
        arr = res_arr[:, tag]
        for i in range(1, cpus):
            arr = numpy.column_stack((arr, res_arr[:, i*len(cutoff_scalar)+tag]))

        e1mc = tool_box.data_fit(fg1, arr[0], arr[1])
        mc1.append(e1mc)
        e2mc = tool_box.data_fit(fg2, arr[3], arr[4])
        mc2.append(e2mc)

        mc = numpy.array([e1mc, e2mc])
        data_path = result_path + "cuts/sym/%s/"%filter_name + cut + "/" + str(round(cut_s,4))+".npz"
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

        g1_extr = [fg1.min(), fg1.max()]
        g2_extr = [fg2.min(), fg2.max()]
        g1_plt_s, g1_plt_e = fg1.min() - (fg1.max() - fg1.min())*0.2, fg1.max() + (fg1.max() - fg1.min())*0.2
        g2_plt_s, g2_plt_e = fg2.min() - (fg2.max() - fg2.min()) * 0.2, fg2.max() + (fg2.max() - fg2.min()) * 0.2
        plt_range = [g1_plt_s, g1_plt_e, g2_plt_s, g2_plt_e]
        pic_path = result_path + "cuts/sym/%s/"%filter_name + cut + "/" + str(round(cut_s,4)) + pic_mc + ".pdf"
        tool_box.mcplot(fg1, arr[0:3,:], fg2, arr[3:6,:], e1mc, e2mc, str(round(cut_s,4)), 'max', plt_range,pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = result_path + "cuts/sym/%s/"%filter_name + cut + "/total.npz"
    numpy.savez(mc_path, mc1, mc2)
    # mc1 = numpy.load(mc_path)['arr_0']
    # mc2 = numpy.load(mc_path)['arr_1']

    x1 = 0
    x2 = 1
    x_coord = [i * 1./cuts_num for i in range(cuts_num)]

    fig = plt.figure(figsize=(10, 10))

    ax1 = fig.add_subplot(221)
    ax1.errorbar(x_coord, mc1[0] - 1, mc1[1], c='coral', capsize=3, label='m1')
    ax1.errorbar(x_coord, mc2[0] - 1, mc2[1], c='royalblue', capsize=3, label='m2')
    ax1.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='green')
    ax1.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax1.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax1.set_xlabel("Cutoff")
    ax1.legend()
    ax1.set_ylabel("m")

    ax2 = fig.add_subplot(222)
    ax2.errorbar(x_coord, mc1[2], mc1[3], c='coral', capsize=3, label='c1')
    ax2.errorbar(x_coord, mc2[2], mc2[3], c='royalblue', capsize=3, label='c2')
    ax2.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='green')
    ax2.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax2.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax2.set_xlabel("Cutoff")
    ax2.legend()
    ax2.set_ylabel("c")

    ax3 = fig.add_subplot(223)
    ax3.hist(cut_data[detected], 100)
    ys = ax3.set_ylim()
    xs = ax3.set_xlim()
    for s in cutoff_scalar:
        ax3.plot([s,s],[ys[0],ys[1]],c="grey",linestyle="--")
    ax3.set_ylim(ys[0], ys[1])
    ax3.set_xlim(xs[0], 1.8*cutoff_scalar[-1]-0.5*cutoff_scalar[-2])

    ax4 = fig.add_subplot(224)
    try:
        detect_num = detected.sum()
        ch = numpy.random.choice([i for i in range(detect_num)], 10000, replace=False)
        ax4.scatter(mag_true[detected][ch], cut_data[detected][ch], s=0.1)
        ax4.set_xlabel("TRUE MAG")
        ax4.set_ylabel(cut)
    except:
        print("PLOT FAILED")
    finally:
        namep = result_path + "cuts/sym/%s/"%filter_name + cut + "/total.pdf"
        plt.suptitle(cut)
        plt.savefig(namep,bbox_inches='tight')
        plt.close()


t2 = time.clock()
if rank == 0:
    log = "END  --  %s, %7s, %9s, %6.2f, %5.2f, %12s"%\
          (total_path.split("/")[-2], filter_name, cut, t2-t1, r_thresh, argv[0])
    logger.info(log)


