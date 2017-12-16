# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import os
import time
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import tool_box
from Fourier_Quad import *
from sys import argv
import numpy
import pandas
import lsstetc

wei, snr_s, snr_e, filter_type, wei_pow, method = argv[1:7]

ts = time.time()

pixel_scale = 0.2
stamp_size = 54

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

shear_input = numpy.load(para_path+"shear.npz")['arr_0']
fg1 = numpy.sort(shear_input[0])
fg2 = numpy.sort(shear_input[1])

# the length of the interval
dfg1 = fg1[1]-fg1[0]
dfg2 = fg2[1]-fg2[0]

# where the result data file are placed
path = result_path + "data/"
# where the result figures will be created
pic_path = result_path + "pic/"

# check the final result cache
final_cache_path = path + 'final_cache.npz'
data_cache_path = path + 'data_cache.npz'
exist = os.path.exists(final_cache_path)
print(path,"\n",pic_path,"\n",final_cache_path,"\n", data_cache_path)
fq = Fourier_Quad(56, 12243)
if exist:
    # 0: to use the result cache data existed to plot the line and estimate the bias, 'm' and 'c'
    # 1: run the program to classify the data ( it could be skipped if this result cache exists) and estimate the shear
    print('0: use the result cache data existed\n1: overwrite it')
    comm = int(input('0/1?:'))
else:
    print('no result cache exists')

if not exist or comm == 1:
    # check the classification cache
    data_cache_exist = os.path.exists(data_cache_path)
    if data_cache_exist:
        print('0: use the classification cache data existed\n1: overwrite it')
        data_comm = int(input('0/1?:'))
    else:
        print('no classification cache')

    if not data_cache_exist or data_comm == 1:
        print('Starting>>>>')
        files = os.listdir(path)
        paths = []
        for i in files:
            if ".npz" and 'chip' in i:
                paths.append(path + i)
        # collect the data from the files and put into 'data_list'
        # 'data' is the final array that contains all the data
        tc1 = time.time()
        data_list = tool_box.classify(paths, 10)
        # cache
        data = data_list[0]
        for k in range(1, len(data_list)):
            data = numpy.row_stack((data, data_list[k]))
        numpy.savez(data_cache_path, data)
        tc2 = time.time()
        print("Classification finished within %.3f <<<<"%(tc2-tc1))

    else:
        print("Loading data cache>>>")
        data = numpy.load(data_cache_path)['arr_0']

    # print("Calculate shear")
    # print(data.shape)
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
    print("area: %d ~ %d" %(numpy.min(area), numpy.max(area)))

    # flux
    flux = data[:, 17]/noise_sig
    print("flux: %.2f ~ %.2f "%(numpy.min(flux), numpy.max(flux)))

    # peak
    peak = data[:, 18]/noise_sig
    print("peak: %.2f ~ %.2f "%(numpy.min(peak), numpy.max(peak)))

    # fsnr
    fsnr = data[:, 19]
    print("fsnr: %.2f ~ %.2f "%(numpy.min(fsnr), numpy.max(fsnr)))

    # snr
    snr = data[:, 20]
    print("snr: %.2f ~ %.2f "%(numpy.min(snr), numpy.max(snr)))

    # input g1
    tag1 = data[:, 4]

    # input g2
    tag2 = data[:, 9]

    select = {"peak": peak, "area": area, "flux": flux, "fsnr": fsnr, "snr": snr}

    # the first 4 rows are the ellipticity,
    # the second 4 rows are the corresponding error bar,
    # the third 4 rows are the corresponding number of samples.
    res_arr1 = numpy.zeros((12, len(fg1)))
    res_arr2 = numpy.zeros((12, len(fg2)))

    for i in range(len(fg1)):
        print(fg1[i])
        idx11 = tag1 > fg1[i] - 0.000001
        idx12 = tag1 < fg1[i] + 0.000001
        ssnr = select[wei][idx11 & idx12]
        idxs = ssnr >= snr_cut_s
        idxe = ssnr <= snr_cut_e
        for na in range(3, 4):
            if na != 3:
                ellip1 = data[:, na]
                idx13 = ellip1 != -10
                r_thresh = r_factor[na]
                idx_r1 = r_thresh > 0.333
                # the measured e1
                e1 = ellip1[idx11&idx12&idx13&idxs&idxe&idx_r1]
                response1 = 2 - measured_es[na][idx11&idx12&idx13&idxs&idxe&idx_r1]
                num1 = len(e1)
                if na!=0:
                    g1_h = numpy.mean(e1)/2#/numpy.mean(measured_esq1)
                else:
                    g1_h = numpy.mean(e1)
                g1_h_sig = numpy.std(e1)/numpy.sqrt(num1)

            else:
                G1 = FG1[idx11&idx12][idxs&idxe]
                N1 = FN[idx11&idx12][idxs&idxe]
                U1 = FU[idx11&idx12][idxs&idxe]

                weight1 = ssnr[idxs&idxe]**wei_pow
                if wei_pow == 0:
                    weight1 = 1

                num1 = len(G1)
                if method == 'sym':
                    g1_h, g1_h_sig = fq.fmin_g(G1, N1, U1, mode=1, bin_num=8)
                else:
                    g1_h = numpy.mean(G1 * weight1) / numpy.mean(N1 * weight1)
                    g1_h_sig = numpy.sqrt(numpy.mean((G1 * weight1)**2)/(numpy.mean(N1 * weight1))**2)/numpy.sqrt(num1)

            res_arr1[na, i] = g1_h
            res_arr1[na+4, i] = g1_h_sig
            res_arr1[na+8, i] = num1
            print(num1, "g1: %.4f, deviation: %.4f * e-4, sig: %.6f" %(g1_h, 10000*(g1_h - fg1[i]), g1_h_sig))
    print('\n')
    for i in range(len(fg2)):
        print(fg2[i])
        idx21 = tag2 > fg2[i] - 0.000001
        idx22 = tag2 < fg2[i] + 0.000001
        ssnr = select[wei][idx21 & idx22]
        idxs = ssnr >= snr_cut_s
        idxe = ssnr <= snr_cut_e
        for na in range(3, 4):
            if na != 3:
                ellip2 = data[:, na+5]
                idx23 = ellip2 != -10
                r_thresh = r_factor[na]
                idx_r2 = r_thresh > 0.3
                # the measured e1
                e2 = ellip2[idx21&idx22&idx23&idxs&idxe&idx_r2]
                response2 = 2 - measured_es[na][idx21&idx22&idx23&idxs&idxe&idx_r2]
                num2 = len(e2)
                if na!=0:
                    g2_h = numpy.mean(e2)/2#/numpy.mean(measured_esq1)
                else:
                    g2_h = numpy.mean(e2)
                g2_h_sig = numpy.std(e2)/numpy.sqrt(num2)

            else:
                G2 = FG2[idx21&idx22][idxs&idxe]
                N2 = FN[idx21&idx22][idxs&idxe]
                U2 = FU[idx21&idx22][idxs&idxe]

                weight2 = ssnr[idxs&idxe]**wei_pow
                if wei_pow == 0:
                    weight2 = 1

                num2 = len(G2)
                if method == 'sym':
                    g2_h, g2_h_sig = fq.fmin_g(G2, N2, U2, mode=2, bin_num=8)
                else:
                    g2_h = numpy.mean(G2 * weight2)/numpy.mean(N2 * weight2)
                    g2_h_sig = numpy.sqrt(numpy.mean((G2 * weight2)**2)/(numpy.mean(N2 * weight2))**2)/numpy.sqrt(num2)
            res_arr2[na, i] = g2_h
            res_arr2[na+4, i] = g2_h_sig
            res_arr2[na+8, i] = num2
            print(num2, "g2: %.5f,  deviation: %.4f * e-4, sig: %.6f" % (g2_h, 10000*(g2_h - fg2[i]), g2_h_sig))

    numpy.savez(final_cache_path, res_arr1, res_arr2)

else:
    text = numpy.load(final_cache_path)
    res_arr1 = text['arr_0']
    res_arr2 = text['arr_1']
tm =time.time()

# fit the line

name = ['KSB', 'BJ02', 'Re-Gaussianization', 'Fourier_Quad']

mc_data_path = path + 'mc_data.xlsx'
mc_data = numpy.zeros((32, 1))
for i in range(3, 4):

    e1mc = tool_box.data_fit(fg1, res_arr1[i], res_arr1[i+4])
    e2mc = tool_box.data_fit(fg2, res_arr2[i], res_arr2[i+4])

    print(e1mc)
    print(e2mc)
    fig = plt.figure(figsize=(20, 10))
    ax = fig.add_subplot(121)

    # xmajorLocator = MultipleLocator(0.005)
    # xmajorFormatter = FormatStrFormatter('%1.3f')
    # xminorLocator = MultipleLocator(0.001)
    #
    # ymajorLocator   = MultipleLocator(0.005)
    # ymajorFormatter = FormatStrFormatter('%1.3f')
    # yminorLocator   = MultipleLocator(0.001)
    #
    # ax.xaxis.set_major_locator(xmajorLocator)
    # ax.xaxis.set_major_formatter(xmajorFormatter)
    #
    # ax.yaxis.set_major_locator(ymajorLocator)
    # ax.yaxis.set_major_formatter(ymajorFormatter)
    #
    # ax.xaxis.set_minor_locator(xminorLocator)
    # ax.yaxis.set_minor_locator(yminorLocator)

    ax.errorbar(fg1, res_arr1[i, :], res_arr1[i + 4, :], ecolor='black', elinewidth='1', fmt='none',capsize=2)
    ax.plot(fg1, e1mc[0] * fg1 + e1mc[2], label=name[i], color='red')
    ax.plot([-0.1, 0.1], [-0.1, 0.1], label='y=x', color='blue')
    ax.scatter(fg1, res_arr1[i, :], c='black')
    for j in range(len(fg1)):
        ax.text(fg1[j], res_arr1[i, j], str(round(res_arr1[i + 8, j] / 1000, 1)) + "K", color="red")

    ax.text(0.1, 0.85, 'm=' + str(round(e1mc[0] - 1, 6))+'$\pm$'+str(round(e1mc[1], 6)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.8, 'c=' + str(round(e1mc[2], 6))+'$\pm$'+str(round(e1mc[3], 6)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.75, "[ " + snr_s + ", " + snr_e + "]", color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g1', fontsize=20)
    plt.ylabel('Est  g1', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.04, 0.04)
    plt.xlim(-0.04, 0.04)
    print('plotted g1')

    #plot g2 line
    ax = fig.add_subplot(122)

    # xmajorLocator = MultipleLocator(0.01)
    # xmajorFormatter = FormatStrFormatter('%1.3f')
    # xminorLocator = MultipleLocator(0.002)
    #
    # ymajorLocator   = MultipleLocator(0.01)
    # ymajorFormatter = FormatStrFormatter('%1.3f')
    # yminorLocator   = MultipleLocator(0.002)
    #
    # ax.xaxis.set_major_locator(xmajorLocator)
    # ax.xaxis.set_major_formatter(xmajorFormatter)
    #
    # ax.yaxis.set_major_locator(ymajorLocator)
    # ax.yaxis.set_major_formatter(ymajorFormatter)
    #
    # ax.xaxis.set_minor_locator(xminorLocator)
    # ax.yaxis.set_minor_locator(yminorLocator)

    ax.errorbar(fg2, res_arr2[i, :], res_arr2[i + 4, :], ecolor='black', elinewidth='1', fmt='none',capsize =2)
    ax.plot(fg2, e2mc[0] * fg2 + e2mc[2], label=name[i], color='red')
    ax.plot([-0.1, 0.1], [-0.1, 0.1], label='y=x', color='blue')
    ax.scatter(fg2, res_arr2[i, :], c='black')
    for j in range(len(fg2)):
        ax.text(fg2[j], res_arr2[i, j], str(round(res_arr2[i + 8, j] / 1000, 1)) + "K", color="red")
    ax.text(0.1, 0.85, 'm=' + str(round(e2mc[0] - 1, 6))+'$\pm$'+str(round(e2mc[1], 6)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.8, 'c=' + str(round(e2mc[2], 6))+'$\pm$'+str(round(e2mc[3], 6)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.75, "[ " + snr_s + ", " + snr_e + "]", color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g2', fontsize=20)
    plt.ylabel('Est  g2', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.04, 0.04)
    plt.xlim(-0.04, 0.04)
    nm = pic_path + name[i] + ".png"
    print(nm)

    print('plotted g2')
    plt.savefig(nm)
    # m1, m2
    mc_data[i] = e1mc[0]
    mc_data[i+16] = e2mc[0]
    # delta_m1, delta_m2
    mc_data[i+4] = e1mc[1]
    mc_data[i+20] = e2mc[1]
    # c1, c2
    mc_data[i+8] = e1mc[2]
    mc_data[i+24] = e2mc[2]
    # delta_c1, delta_c2
    mc_data[i+12] = e1mc[3]
    mc_data[i+28] = e2mc[3]

if filter_type != 'none':
    mc_col = [filter_type]
    if os.path.exists(mc_data_path):
        df = pandas.read_excel(mc_data_path)
        df[mc_col[0]] = mc_data
        df.to_excel(mc_data_path)
    else:
        dex = ['Km1', 'Bm1', 'Rm1', 'Fm1', 'Kdm1', 'Bdm1', 'Rdm1', 'Fdm1', 'Kc1', 'Bc1', 'Rc1', 'Fc1', 'Kdc1', 'Bdc1', 'Rdc1', 'Fdc1','Km2', 'Bm2', 'Rm2', 'Fm2', 'Kdm2', 'Bdm2', 'Rdm2', 'Fdm2', 'Kc2', 'Bc2', 'Rc2', 'Fc2', 'Kdc2', 'Bdc2', 'Rdc2', 'Fdc2']
        mc_df = pandas.DataFrame(data=mc_data, index=dex, columns=mc_col)
        mc_df.to_excel(mc_data_path)
te = time.time()

print(tm-ts, te-tm)

