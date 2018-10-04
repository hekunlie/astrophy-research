# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import time
import tool_box
from Fourier_Quad import *
from sys import argv
import numpy
from mpi4py import MPI
import lsstetc
import h5py


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.clock()

ch, snr_s, snr_e, scale, del_bin = argv[1:8]

pixel_scale = 0.2
stamp_size = 60

snr_cut_s = float(snr_s)
snr_cut_e = float(snr_e)
del_bin = int(del_bin)

ini_path = "%s/work/envs/envs.dat"%my_home
path_items = tool_box.config(ini_path,['get','get','get','get'], [['selection_bias', "ptsm_path", '1'],
                                                                  ['selection_bias', "ptsm_path_result", '1'],
                                                                  ['selection_bias', "ptsm_path_para", '1'],
                                                                  ['selection_bias', "ptsm_path_pic", '1']])

total_path = path_items[0]
result_path = path_items[1]
para_path = path_items[2]
pic_path = path_items[3]

shear_input = numpy.load(para_path+"shear.npz")
fg1 = shear_input['arr_0']
fg2 = shear_input['arr_1']

# where the result data file are placed
path = result_path + "data/"

final_cache_path = path + '%d_%s_%s_final_cache.npz'%(del_bin, ch, snr_s)
for s in range(int(scale)):
    if scale == 1:
        data_cache_path = path + "data_%d.hdf5"%rank
    else:
        data_cache_path = path + 'data_%d_%d.hdf5'%(rank, s)
    f = h5py.File(data_cache_path, 'r')
    if s == 0:
        data = f["/data"].value.copy()
    else:
        data = numpy.row_stack((data, f["/data"].value))
    f.close()

fq = Fourier_Quad(stamp_size, 123)

flux = data[:, 7]
detect = flux > 0

# F_Q data
FG1 = data[:, 2]
FG2 = data[:, 3]
FN = data[:, 4]
FU = data[:, 5]
FV = data[:, 6]

prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
noise_sig = prop.sigma_sky

# flux
flux = data[:, 8]/noise_sig
# half_light_flux
hflux = data[:, 8]/noise_sig
# peak
peak = data[:, 9]/noise_sig
# area
area = data[:, 12]
# half_light_area
harea = data[:, 11]
# snr
snr = data[:, 7]
# flux2
flux2 = data[:, 10]
# flux_alt
flux_alt = data[:, 9]


if rank == 0:
    print("flux: %.2f ~ %.2f\n"%(numpy.min(flux), numpy.max(flux)))
    print("peak: %.2f ~ %.2f\n"%(numpy.min(peak), numpy.max(peak)))
    print("flux2: %.2f ~ %.2f\n"%(numpy.min(flux2), numpy.max(flux2)))
    print("osnr: %.2f ~ %.2f\n"%(numpy.min(snr), numpy.max(snr)))
    print("area: %d ~ %d\n" % (numpy.min(area), numpy.max(area)))

# sex_path = total_path + "result/data/sex25_%d_1.5.npz"%rank
# sex_data = numpy.load(sex_path)["arr_0"]
# sex_snr = sex_data[:, 0]
# sex_path = path + "sex25_%d_1.5.npz"%rank
# sex_snr = numpy.load(sex_path)['arr_0'][:, 0]

select = {"peak": peak, "flux": flux, "flux2": flux2, "area": area, "snr": snr,
          "hflux": hflux, "harea": harea, "flux_alt": flux_alt}#, "sex_snr": sex_snr}

res_arr = numpy.zeros((3, 2))
sp = res_arr.shape

ssnr = select[ch]
idxs = ssnr > snr_cut_s
idxe = ssnr <= snr_cut_e


G1 = FG1[idxs&idxe]
G2 = FG2[idxs & idxe]
N = FN[idxs&idxe]
U = FU[idxs&idxe]
DE1 = N + U
DE2 = N - U

num = len(G1)
print(rank, num)

g1_xi2_pic = pic_path + "%d_%s_%d_g1_xi2.png"%(rank, ch, del_bin)
g2_xi2_pic = pic_path + "%d_%s_%d_g2_xi2.png"%(rank, ch, del_bin)
g1_h, g1_h_sig = fq.fmin_g_new(G1, DE1, bin_num=8, ig_num=del_bin, pic_path=g1_xi2_pic)
g2_h, g2_h_sig = fq.fmin_g_new(G2, DE2, bin_num=8, ig_num=del_bin, pic_path=g2_xi2_pic)

# else:
#     weight = 1
#     g1_h = numpy.mean(G1 * weight) / numpy.mean(N * weight)
#     g1_h_sig = numpy.sqrt(numpy.mean((G1 * weight)**2)/(numpy.mean(N * weight))**2)/numpy.sqrt(num)
#
#     g2_h = numpy.mean(G2 * weight) / numpy.mean(N * weight)
#     g2_h_sig = numpy.sqrt(numpy.mean((G2 * weight) ** 2) / (numpy.mean(N * weight)) ** 2) / numpy.sqrt(num)

res_arr[0] = g1_h, g2_h
res_arr[1] = g1_h_sig, g2_h_sig
res_arr[2] = num, num

res_arr = [g1_h, g1_h_sig, num, g2_h, g2_h_sig, num]
res_arr = comm.gather(res_arr, root=0)

# if rank > 0:
#     comm.Send(res_arr, dest=0, tag=rank)
# else:
#     for procs in range(1, cpus):
#         recvs = numpy.empty(sp, dtype=numpy.float64)
#         comm.Recv(recvs, source=procs, tag=procs)
#         res_arr = numpy.column_stack((res_arr, recvs))
if rank == 0:
    res_arr = numpy.array(res_arr)
    res_arr = res_arr.T
    # fit the line
    for i in range(3, 4):
        arr1 = res_arr[:3]
        arr2 = res_arr[3:]
        e1mc = tool_box.data_fit(fg1, arr1[0], arr1[1])
        e2mc = tool_box.data_fit(fg2, arr2[0], arr2[1])
        numpy.savez(final_cache_path,arr1, arr2, e1mc, e2mc)
        for p in range(cpus):
            print("num: %4.1f W, g1: %8.5f, m_g1: %8.5f, sig: %8.5f, devi: %4.2f * e^-4, shape noise: %6.4f"
                  %(arr1[2,p]/10000, fg1[p], arr1[0,p], arr1[1,p], 10000*(arr1[0,p]-fg1[p]), numpy.sqrt(arr1[2,p])*arr1[1,p]))
        print("\n")
        for p in range(cpus):
            print("num: %4.1f W, g2: %8.5f, m_g2: %8.5f, sig: %8.5f, devi: %4.2f * e^-4, shape noise: %6.4f "
                  %(arr2[2,p]/10000, fg2[p], arr2[0,p], arr2[1,p], 10000*(arr2[0,p]-fg2[p]), numpy.sqrt(arr2[2,p])*arr2[1,p]))

        if e1mc[0]-1 - 2*e1mc[1] < 0 < e1mc[0]-1 + 2*e1mc[1]:
            m1_b = "no m1 bias"
        else:
            m1_b = "m1 bias"
        if e2mc[0]-1 - 2*e2mc[1] < 0 < e2mc[0]-1 + 2*e2mc[1]:
            m2_b = "no m2 bias"
        else:
            m2_b = "m2 bias"

        if e1mc[2] - 2*e1mc[3] < 0 < e1mc[2] + 2*e1mc[3]:
            c1_b = "no c1 bias"
        else:
            c1_b = "c1 bias"
        if e2mc[2] - 2*e2mc[3] < 0 < e2mc[2] + 2*e2mc[3]:
            c2_b = "no c2 bias"
        else:
            c2_b = "c2 bias"
        print("%10s: %8.5f (%6.5f), %10s: %10.6f (%.6f)"%(m1_b, e1mc[0]-1, e1mc[1], c1_b, e1mc[2], e1mc[3]))
        print("%10s: %8.5f (%6.5f), %10s: %10.6f (%.6f)"%(m2_b, e2mc[0]-1, e2mc[1], c2_b, e2mc[2], e2mc[3]))
        nm = pic_path + "mc_%s.png"%ch
        tool_box.mcplot(fg1, arr1, fg2, arr2, e1mc, e2mc, snr_s, 'max',[-0.03,0.03,-0.03,0.03], nm)

te = time.clock()
if rank == 0:
    print(te-ts)
