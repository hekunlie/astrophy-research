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
# rfactor_path = result_path + "data/resolution_factor_%d.npz"%rank

if rank == 0:
    log = "START -- %s, %7s, %9s,      , %5.2f, %12s"%(total_path.split("/")[-2], filter_name, cut, r_thresh, argv[0])
    logger = tool_box.get_logger("%s/work/selection_bias/sym_mc_plot/cutoff.dat"%my_home)
    logger.info(log)

shear_esti_h5path = result_path + "data/data_%d.hdf5"%rank
shear_esti_file = h5py.File(shear_esti_h5path, "r")
shear_esti_data = shear_esti_file["/data"].value
shear_esti_file.close()

fq_snr_h5path = result_path + "data/data_%s/data_%d.hdf5" % (sig, rank)
fq_snr_file = h5py.File(fq_snr_h5path, "r")
set_name = list(fq_snr_file.keys())[0]
fq_snr_data = fq_snr_file[set_name].value
fq_snr_file.close()

fq = Fourier_Quad(60, 152356)
noise_sig = 65

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

detect_ = fq_snr_data[:, 0]
detected = detect_ > 0

cuts_num = 10

# snr
snr = fq_snr_data[:, 0]
snr_sort = numpy.sort(snr[detected])
snr_step = int(len(snr_sort)/cuts_num)
snr_cut = [snr_sort[i*snr_step] for i in range(cuts_num)]

# flux
flux = fq_snr_data[:, 1]/noise_sig
flux_sort = numpy.sort(flux[detected])
flux_step = int(len(flux_sort)/cuts_num)
flux_cut = [flux_sort[i*flux_step] for i in range(cuts_num)]

# flux_alt
flux_alt = fq_snr_data[:, 2]
flux_alt_sort = numpy.sort(flux_alt[detected])
flux_alt_step = int(len(flux_alt_sort)/cuts_num)
flux_alt_cut = [flux_alt_sort[i*flux_alt_step] for i in range(cuts_num)]

# flux2
flux2 = fq_snr_data[:, 3]
flux2_sort = numpy.sort(flux2[detected])
flux2_step = int(len(flux2_sort)/cuts_num)
flux2_cut = [flux2_sort[i*flux2_step] for i in range(cuts_num)]

# great3 snr
gsnr = fq_snr_data[:, 4]
gsnr_sort = numpy.sort(gsnr[detected])
gsnr_step = int(len(gsnr_sort)/cuts_num)
gsnr_cut = [gsnr_sort[i*gsnr_step] for i in range(cuts_num)]

# area
area = fq_snr_data[:, 5]
area_sort = numpy.sort(area[detected])
area_step = int(len(area_sort)/cuts_num)
area_cut = [area_sort[i*area_step] for i in range(cuts_num)]


sex_data = numpy.load(sex_path)["arr_0"]

sex_snr = sex_data[:, 0]
sex_idx = sex_snr > 0
sex_snr_sort = numpy.sort(sex_snr[sex_idx])
sex_snr_step = int(len(sex_snr_sort)/cuts_num)
sex_snr_cut = [sex_snr_sort[i*sex_snr_step] for i in range(cuts_num)]

flux_auto = sex_data[:, 1]
flux_err = sex_data[:, 2]
err_idx = flux_err == 0
sex_undetect = sex_snr == 0
if rank == 0:
    print(sig,err_idx.sum()-sex_undetect.sum())
flux_err[err_idx] = 1
snr_auto = flux_auto/flux_err
snr_auto_sort = numpy.sort(snr_auto[sex_idx])
snr_auto_step = int(len(snr_auto_sort)/cuts_num)
snr_auto_cut = [snr_auto_sort[i*snr_auto_step] for i in range(cuts_num)]

mag_auto = sex_data[:, 3]
mag_auto_sort = numpy.sort(-mag_auto[sex_idx])
mag_auto_step = int(len(mag_auto_sort)/cuts_num)
mag_auto_cut = [mag_auto_sort[i*mag_auto_step] for i in range(cuts_num)]

sex_area = sex_data[:, 4]
sex_area_sort = numpy.sort(sex_area[sex_idx])
sex_area_step = int(len(sex_area_sort)/cuts_num)
sex_area_cut = [sex_area_sort[i*sex_area_step] for i in range(cuts_num)]

# # Resolution factor
# rfactor = numpy.load(rfactor_path)["arr_0"][:,0] # [:,0] to avoid memory error
# idx_rf = rfactor >= r_thresh

select = {"snr":     (snr, snr_cut),          "flux":     (flux, flux_cut),
          "area":    (area, area_cut),        "flux2":   (flux2, flux2_cut),
          "flux_alt": (flux_alt, flux_alt_cut),
          "sex_snr": (sex_snr, sex_snr_cut),  "sex_area": (sex_area, sex_area_cut),
          "snr_auto": (snr_auto, snr_auto_cut),  "mag_auto": (-mag_auto, mag_auto_cut)}

res_arr = numpy.zeros((6, len(select[cut][1])))


detected_label = {"flux": detected, "hflux": detected, "peak": detected, "area": detected, "harea": detected,
                  "snr": detected, "flux2": detected, "flux_alt": detected, "sex_snr": sex_idx, "sex_area": sex_idx,
                  "snr_auto": sex_idx, "mag_auto": sex_idx}

for tag, cut_s in enumerate(select[cut][1]):
    idx = select[cut][0] >= cut_s
    num = len(MG1[detected_label[cut]&idx])

    nm1 = MG1[detected_label[cut]&idx]
    de1 = DE1[detected_label[cut]&idx]
    nm2 = MG2[detected_label[cut]&idx]
    de2 = DE2[detected_label[cut]&idx]

    g1_h, g1_sig = fq.fmin_g_new(nm1, de1, bin_num=8)
    g2_h, g2_sig = fq.fmin_g_new(nm2, de2, bin_num=8)

    # weight = select[cut][0][idx]
    # g1_h = numpy.mean(mg1[idx] / weight) / numpy.mean(mn[idx] / weight)
    # g1_sig = numpy.sqrt(numpy.mean((mg1[idx] / weight) ** 2) / (numpy.mean(mn[idx] / weight)) ** 2) / numpy.sqrt(num)
    # g2_h = numpy.mean(mg2[idx] / weight) / numpy.mean(mn[idx] / weight)
    # g2_sig = numpy.sqrt(numpy.mean((mg2[idx] / weight) ** 2) / (numpy.mean(mn[idx] / weight)) ** 2) / numpy.sqrt(num)

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
        data_path = total_path + "result/cuts/sym/%s/"%filter_name + cut + "/" + str(round(cut_s,4))+".npz"
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

        pic_path = total_path + "result/cuts/sym/%s/"%filter_name + cut + "/" + str(round(cut_s,4)) + pic_mc + ".eps"
        tool_box.mcplot(fg1, arr[0:3,:], fg2, arr[3:6,:], e1mc, e2mc, str(round(cut_s,4)), 'max', [-0.03,0.03,-0.03,0.03],pic_path)
        pic_path = total_path + "result/cuts/sym/%s/"%filter_name + cut + "/" + str(round(cut_s,4)) + pic_mc + ".png"
        tool_box.mcplot(fg1, arr[0:3, :], fg2, arr[3:6, :], e1mc, e2mc, str(round(cut_s,4)), 'max',[-0.03,0.03,-0.03,0.03],pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = total_path + "result/cuts/sym/%s/"%filter_name + cut + "/total.npz"
    numpy.savez(mc_path, mc1, mc2)
    # mc1 = numpy.load(mc_path)['arr_0']
    # mc2 = numpy.load(mc_path)['arr_1']

    x1 = 0
    x2 = 1
    x_coord = [i * 1./cuts_num for i in range(cuts_num)]

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

    namep = total_path + "result/cuts/sym/%s/"%filter_name + cut + "/total.eps"
    plt.suptitle(cut)
    plt.savefig(namep)
    plt.close()

t2 = time.clock()
if rank == 0:
    log = "END  --  %s, %7s, %9s, %6.2f, %5.2f, %12s"%\
          (total_path.split("/")[-2], filter_name, cut, t2-t1, r_thresh, argv[0])
    logger.info(log)


