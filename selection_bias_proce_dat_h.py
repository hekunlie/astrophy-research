# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from sys import path
path.append('/home/hkli/work/fourier_quad')
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

wei, snr_s, snr_e, wei_pow, method, scale, del_bin = argv[1:8]

pixel_scale = 0.2
stamp_size = 90

snr_cut_s = float(snr_s)
snr_cut_e = float(snr_e)
wei_pow = int(wei_pow)
del_bin = int(del_bin)

with open("/home/hkli/work/envs/envs.dat", "r") as f:
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
pic_path = result_path + "pic/%d_%s_"%(del_bin,wei)

final_cache_path = path + '%d_%s_%s_final_cache.npz'%(del_bin, wei,snr_s)
for s in range(int(scale)):
    if scale == 1:
        data_cache_path = path + "data_%d.hdf5"%rank
    else:
        data_cache_path = path + 'data_f_%d_%d.hdf5'%(rank,s)
    f = h5py.File(data_cache_path,'r')
    if s == 0:
        data = f["/data"].value
    else:
        data = numpy.row_stack((data, f["/data"].value))
    f.close()
sex_path = path + "sex25_%d_1.5.npz"%rank
sex_snr = numpy.load(sex_path)['arr_0'][:, 0]

fq = Fourier_Quad(stamp_size, 123)

# F_Q data
FG1 = data[:, 2]
FG2 = data[:, 3]
FN = data[:, 4]
FU = data[:, 5]
FV = data[:, 6]

prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
noise_sig = prop.sigma_sky

# osnr
osnr = data[:, 7]
# flux
flux = data[:, 8]/noise_sig
# peak
peak = data[:, 9]/noise_sig
# fsnr
fsnr = data[:, 10]
# snr
snr = data[:, 11]
# area
area = data[:, 12]

if rank == 0:
    print("flux: %.2f ~ %.2f\n"%(numpy.min(flux), numpy.max(flux)))
    print("peak: %.2f ~ %.2f\n"%(numpy.min(peak), numpy.max(peak)))
    print("snr: %.2f ~ %.2f\n"%(numpy.min(snr), numpy.max(snr)))
    print("fsnr: %.2f ~ %.2f\n"%(numpy.min(fsnr), numpy.max(fsnr)))
    print("osnr: %.2f ~ %.2f\n"%(numpy.min(osnr), numpy.max(osnr)))
    print("area: %d ~ %d\n" % (numpy.min(area), numpy.max(area)))

select = {"peak": peak, "flux": flux, "fsnr": fsnr, "snr": snr, "area": area, "sex": sex_snr}

res_arr = numpy.zeros((3, 2))
sp = res_arr.shape

ssnr = select[wei]
idxs = ssnr >= snr_cut_s
idxe = ssnr <= snr_cut_e


G1 = FG1[idxs&idxe]
G2 = FG2[idxs & idxe]
N = FN[idxs&idxe]
U = FU[idxs&idxe]
weight = ssnr[idxs&idxe]**wei_pow
if wei_pow == 0:
    weight = 1

num = len(G1)
if method == 'sym':
    g1_xi2_pic = pic_path + "%d_g1_xi2.png"%rank
    g2_xi2_pic = pic_path + "%d_g2_xi2.png"%rank
    g1_h, g1_h_sig = fq.fmin_g(G1, N, U, mode=1, bin_num=18, ig_num=del_bin, pic_path=g1_xi2_pic)
    g2_h, g2_h_sig = fq.fmin_g(G2, N, U, mode=2, bin_num=18, ig_num=del_bin, pic_path=g2_xi2_pic)
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
        print(int(del_bin),"%10s: %8.5f (%6.5f), %10s: %10.6f (%.6f)"%(m1_b, e1mc[0]-1, e1mc[1], c1_b, e1mc[2], e1mc[3]))
        print(int(del_bin),"%10s: %8.5f (%6.5f), %10s: %10.6f (%.6f)"%(m2_b, e2mc[0]-1, e2mc[1], c2_b, e2mc[2], e2mc[3]))
        nm = pic_path + name[i] + ".png"
        tool_box.mcplot(fg1, arr1, fg2, arr2, e1mc, e2mc, snr_s, 'max', nm)

te = time.clock()
if rank == 0:
    print(te-ts)