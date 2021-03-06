import matplotlib
matplotlib.use("Agg")
import numpy
import matplotlib.pyplot as plt
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
import datetime


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

cut, file_num = argv[1], int(argv[2])

t1 = time.clock()
source = "ptsm"
env_path = "%s/work/envs/envs.dat"%my_home
path_items = tool_box.config(env_path, ['get','get','get'], [['selection_bias', "%s_path"%source, '1'],
                                                             ['selection_bias', "%s_path_result"%source, '1'],
                                                             ['selection_bias', "%s_path_para"%source, '1']])

total_path, result_path, para_path = path_items
shear_path = para_path + "shear.npz"
shear = numpy.load(shear_path)
fg1 = shear["arr_0"]
fg2 = shear["arr_1"]
sex_path = total_path + "result/data/sex2_1.5_%d.npz"%rank
for i in range(file_num):
    h5path = result_path + "data/data_%d_%d.hdf5"%(rank,i)
    f = h5py.File(h5path, "r")
    temp = f["/data"].value
    if i == 0:
        data = temp.copy()
    else:
        data = numpy.row_stack((data, temp))
    f.close()

fq = Fourier_Quad(60, 152356)
noise_sig = 380.64

detect_ = data[:, 7]
detected = detect_ > 0

cuts_num = 20

# flux
flux = data[:, 7]/noise_sig
flux_sort = numpy.sort(flux[detected])
flux_step = int(len(flux_sort)/cuts_num)
flux_cut = [flux_sort[i*flux_step] for i in range(cuts_num)]

# half_light_flux
hflux = data[:, 8]/noise_sig
hflux_sort = numpy.sort(hflux[detected])
hflux_step = int(len(hflux_sort)/cuts_num)
hflux_cut = [hflux_sort[i*hflux_step] for i in range(cuts_num)]

# peak
peak = data[:, 9]/noise_sig
peak_sort = numpy.sort(peak[detected])
peak_step = int(len(peak_sort)/cuts_num)
peak_cut = [peak_sort[i*peak_step] for i in range(cuts_num)]

# area
area = data[:, 10]
area_sort = numpy.sort(area[detected])
area_step = int(len(area_sort)/cuts_num)
area_cut = [area_sort[i*area_step] for i in range(cuts_num)]

# half_light_area
harea = data[:, 11]
harea_sort = numpy.sort(harea[detected])
harea_step = int(len(harea_sort)/cuts_num)
harea_cut = [harea_sort[i*harea_step] for i in range(cuts_num)]

# snr
snr = data[:, 12]
snr_sort = numpy.sort(snr[detected])
snr_step = int(len(snr_sort)/cuts_num)
snr_cut = [snr_sort[i*snr_step] for i in range(cuts_num)]

# flux2
flux2 = data[:, 10]
flux2_sort = numpy.sort(flux2[detected])
flux2_step = int(len(flux2_sort)/cuts_num)
flux2_cut = [flux2_sort[i*flux2_step] for i in range(cuts_num)]

# flux_alt
flux_alt = data[:, 14]
flux_alt_sort = numpy.sort(flux_alt[detected])
flux_alt_step = int(len(flux_alt_sort)/cuts_num)
flux_alt_cut = [flux_alt_sort[i*flux_alt_step] for i in range(cuts_num)]



sex_data = numpy.load(sex_path)["arr_0"]

sex_snr = sex_data[:, 0]
sex_idx = sex_snr > 0
sex_snr_sort = numpy.sort(sex_snr[sex_idx])
sex_snr_step = int(len(sex_snr_sort)/cuts_num)
sex_snr_cut = [sex_snr_sort[i*sex_snr_step] for i in range(cuts_num)]

sex_area = sex_data[:, 1]
sex_area_sort = numpy.sort(sex_area[sex_idx])
sex_area_step = int(len(sex_area_sort)/cuts_num)
sex_area_cut = [sex_area_sort[i*sex_area_step] for i in range(cuts_num)]

mag_iso = sex_data[:, 2]
mag_iso_sort = numpy.sort(-mag_iso[sex_idx])
mag_iso_step = int(len(mag_iso_sort)/cuts_num)
mag_iso_cut = [mag_iso_sort[i*mag_iso_step] for i in range(cuts_num)]

mag_auto = sex_data[:, 3]
mag_auto_sort = numpy.sort(-mag_auto[sex_idx])
mag_auto_step = int(len(mag_auto_sort)/cuts_num)
mag_auto_cut = [mag_auto_sort[i*mag_auto_step] for i in range(cuts_num)]

mag_petro = sex_data[:, 4]
mag_petro_sort = numpy.sort(-mag_petro[sex_idx])
mag_petro_step = int(len(mag_petro_sort)/cuts_num)
mag_petro_cut = [mag_petro_sort[i*mag_petro_step] for i in range(cuts_num)]

mag_win = sex_data[:, 5]
mag_win_sort = numpy.sort(-mag_win[sex_idx])
mag_win_step = int(len(mag_win_sort)/cuts_num)
mag_win_cut = [mag_win_sort[i*mag_win_step] for i in range(cuts_num)]

select = {"snr":     (snr, snr_cut),          "flux":     (flux, flux_cut),
          "hflux":   (hflux, hflux_cut),      "peak":     (peak, peak_cut),
          "area":    (area, area_cut),        "harea":    (harea, harea_cut),
          "flux2":   (flux2, flux2_cut),      "flux_alt": (flux_alt, flux_alt_cut),
          "sex_snr": (sex_snr, sex_snr_cut),  "sex_area": (sex_area, sex_area_cut),
          "mag_iso": (-mag_iso, mag_iso_cut),  "mag_auto": (-mag_auto, mag_auto_cut),
          "mag_petro": (-mag_petro, mag_petro_cut), "mag_win": (-mag_win, mag_win_cut)}

res_arr = numpy.zeros((6, len(select[cut][1])))

MG1 = data[:, 2]
MG2 = data[:, 3]
MN = data[:, 4]
MU = data[:, 5]
MV = data[:, 6]
# be careful that the "MU" defined in FRESH is the different from that in ours
# MN + (-) MU for our definition (g1(2)) of MU and MV which is the same as
# those in the paper Zhang et al. 2017 ApJ, 834:8
DE1 = MN + MU
DE2 = MN - MU


detected_label = {"flux": detected, "hflux": detected, "peak": detected, "area": detected, "harea": detected,
                  "snr": detected, "flux2": detected, "flux_alt": detected, "sex_snr": sex_idx, "sex_area": sex_idx,
                  "mag_iso": sex_idx, "mag_auto": sex_idx, "mag_petro": sex_idx, "mag_win": sex_idx}

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
        data_path = total_path + "result/cuts/sym/" + cut + "/" + str(round(cut_s,4))+".npz"
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

        pic_path = total_path + "result/cuts/sym/" + cut + "/" + str(round(cut_s,4)) + pic_mc + ".eps"
        tool_box.mcplot(fg1, arr[0:3,:], fg2, arr[3:6,:], e1mc, e2mc, str(round(cut_s,4)), 'max', [-0.03,0.03,-0.03,0.03],pic_path)
        pic_path = total_path + "result/cuts/sym/" + cut + "/" + str(round(cut_s,4)) + pic_mc + ".png"
        tool_box.mcplot(fg1, arr[0:3, :], fg2, arr[3:6, :], e1mc, e2mc, str(round(cut_s,4)), 'max',[-0.03,0.03,-0.03,0.03],pic_path)

    mc1 = numpy.array(mc1).T
    mc2 = numpy.array(mc2).T
    mc_path = total_path + "result/cuts/sym/" + cut + "/total.npz"
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

    namep = total_path + "result/cuts/sym/" + cut + "/total.eps"
    plt.savefig(namep)
    plt.close()

t2 = time.clock()
if rank == 0:
    nowTime = datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    print(total_path.split("/")[-1], nowTime, cut, t2-t1)

