# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
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

wei, snr_s, snr_e, wei_pow, method = argv[1:7]

pixel_scale = 0.2
stamp_size = 64

snr_cut_s = int(snr_s)
snr_cut_e = int(snr_e)
wei_pow = int(wei_pow)

with open("/home/hklee/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "total_data" in path:
        total_path = path.split("=")[1]
    elif "result" in path:
        result_path = path.split("=")[1]
    elif "parameter" in path:
        para_path = path.split("=")[1]

shear_input = numpy.load(para_path+"shear.npz")
fg1 = shear_input['arr_0']
fg2 = shear_input['arr_1']

# where the result data file are placed
path = result_path + "data/"
# where the result figures will be created
pic_path = result_path + "pic/"

# check the final result cache
final_cache_path = path + 'final_cache.npz'
data_cache_path = path + 'data_%d.hdf5'%rank

f = h5py.File(data_cache_path,'r')
data = f["/data"].value
f.close()
fq = Fourier_Quad(stamp_size, 52232345)
# correction
KSB_r = data[:, 13]
BJ_r = data[:, 14]
RG_r = data[:, 15]
r_factor = [KSB_r, BJ_r, RG_r]
k_me = data[:, 0]**2 + data[:, 5]**2
b_me = data[:, 1]**2 + data[:, 6]**2
r_me = data[:, 2]**2 + data[:, 7]**2
measured_es = [k_me, b_me, r_me]

# F_Q data
FG1 = data[:, 3]
FG2 = data[:, 8]
FN = data[:, 10]
FU = data[:, 11]
FV = data[:, 12]

prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
noise_sig = prop.sigma_sky

# area
area = data[:, 16]
# flux
flux = data[:, 17]/noise_sig
# peak
peak = data[:, 18]/noise_sig
# snr
snr = data[:, 19]
# fsnr
fsnr = data[:, 20]
# osnr
osnr = data[:, 21]

if rank == 0:
    print("area: %d ~ %d\n" %(numpy.min(area), numpy.max(area)))
    print("flux: %.2f ~ %.2f\n"%(numpy.min(flux), numpy.max(flux)))
    print("peak: %.2f ~ %.2f\n"%(numpy.min(peak), numpy.max(peak)))
    print("snr: %.2f ~ %.2f\n"%(numpy.min(snr), numpy.max(snr)))
    print("fsnr: %.2f ~ %.2f\n"%(numpy.min(fsnr), numpy.max(fsnr)))
    print("osnr: %.2f ~ %.2f\n"%(numpy.min(osnr), numpy.max(osnr)))

# input g1
tag1 = data[:, 4]

# input g2
tag2 = data[:, 9]

select = {"peak": peak, "area": area, "flux": flux, "fsnr": fsnr, "snr": snr}

res_arr = numpy.zeros((3, 2))
sp = res_arr.shape

ssnr = select[wei]
idxs = ssnr >= snr_cut_s
idxe = ssnr <= snr_cut_e

for na in range(3, 4):
    if na != 3:
        ellip1 = data[:, na]
        idx13 = ellip1 != -10
        r_thresh = r_factor[na]
        idx_r1 = r_thresh > 0.333
        # the measured e1
        e1 = ellip1[idx13&idxs&idxe&idx_r1]
        response1 = 2 - measured_es[na][idx13&idxs&idxe&idx_r1]
        num1 = len(e1)
        if na!=0:
            g1_h = numpy.mean(e1)/2#/numpy.mean(measured_esq1)
        else:
            g1_h = numpy.mean(e1)
        g1_h_sig = numpy.std(e1)/numpy.sqrt(num1)

    else:
        G1 = FG1[idxs&idxe]
        G2 = FG2[idxs & idxe]
        N = FN[idxs&idxe]
        U = FU[idxs&idxe]
        weight = ssnr[idxs&idxe]**wei_pow
        if wei_pow == 0:
            weight = 1

        num = len(G1)
        if method == 'sym':
            g1_h, g1_h_sig = fq.fmin_g(G1, N, U, mode=1, bin_num=8)
            g2_h, g2_h_sig = fq.fmin_g(G2, N, U, mode=2, bin_num=8)
        else:
            g1_h = numpy.mean(G1 * weight) / numpy.mean(N * weight)
            g1_h_sig = numpy.sqrt(numpy.mean((G1 * weight)**2)/(numpy.mean(N * weight))**2)/numpy.sqrt(num)

            g2_h = numpy.mean(G2 * weight) / numpy.mean(N * weight)
            g2_h_sig = numpy.sqrt(numpy.mean((G2 * weight) ** 2) / (numpy.mean(N * weight)) ** 2) / numpy.sqrt(num)

    res_arr[0] = g1_h, g2_h
    res_arr[1] = g1_h_sig, g2_h_sig
    res_arr[2] = num, num

    if rank > 0:
        comm.Send(res_arr, dest=0, tag=rank)
    else:
        for procs in range(1, cpus):
            recvs = numpy.empty(sp, dtype=numpy.float64)
            comm.Recv(recvs, source=procs, tag=procs)
            res_arr = numpy.column_stack((res_arr, recvs))

        # fit the line
        name = ['KSB', 'BJ02', 'Re-Gaussianization', "Fourier_Quad"]

        for i in range(3, 4):
            idx = [k for k in range(0, 27, 2)]
            arr1 = res_arr[:, idx]
            idx = [k for k in range(1, 28, 2)]
            arr2 = res_arr[:, idx]
            e1mc = tool_box.data_fit(fg1, arr1[0], arr1[1])
            e2mc = tool_box.data_fit(fg2, arr2[0], arr2[1])
            numpy.savez(path+"final_cache.npz",arr1, arr2, e1mc, e2mc)
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
            nm = pic_path + name[i] + ".png"
            tool_box.mcplot(fg1, arr1, fg2, arr2, e1mc, e2mc, snr_s, 'max', nm)

te = time.clock()
if rank == 0:
    print(te-ts)

