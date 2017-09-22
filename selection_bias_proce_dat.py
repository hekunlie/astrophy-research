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

snr_s, snr_e, filter_type = argv[1:4]

ts = time.time()

snr_cut_s = int(snr_s)
snr_cut_e = int(snr_e)

shear_input = numpy.load('/lmc/selection_bias/shear.npz')
fg1 = shear_input['arr_0']
fg2 = shear_input['arr_1']

# the length of the interval
dfg1 = fg1[1]-fg1[0]
dfg2 = fg2[1]-fg2[0]

# where the result data file are placed
path = '/lmc/selection_bias/result/data/'
# where the result figures will be created
pic_path = '/lmc/selection_bias/result/pic/'

# check the final result cache
exist = os.path.exists(path+'final_cache.npz')
data_cache_path = path + 'data_cache.npz'

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
            if ".xlsx" and 'chip' in i:
                paths.append(path + i)
        # collect the data from the files and put into 'data_list'
        # 'data' is the final array that contains all the data
        tc1 = time.time()
        data_list = tool_box.classify(paths, 16)[0]
        # cache
        # data_cache = shelve.open(path+'data_cache')
        # data_cache["data"] = data_list
        # data_cache.close()
        data = data_list[0]
        for k in range(1, len(data_list)):
            data = numpy.row_stack((data, data_list[k]))
        numpy.savez(data_cache_path, data)
        tc2 = time.time()
        print("Classification finished within %.3f <<<<"%(tc2-tc1))

    else:
        print("Loading data cache>>>")
        data = numpy.load(data_cache_path)['arr_0']

    print("Calculate shear")
    print(data.shape)
    # correction
    KSB_r = data[:, 13]
    BJ_r = data[:, 14]
    RG_r = data[:, 15]
    r_factor = [KSB_r, BJ_r, RG_r]
    k_me = data[:, 0]**2 + data[:, 5]**2
    b_me = data[:, 1]**2 + data[:, 6]**2
    r_me = data[:, 2]**2 + data[:, 7]**2
    measured_es = [k_me, b_me, r_me]

    # flux
    flux = data[:, 17]

    # peak
    peak = data[:, 18]

    # input g1
    tag1 = data[:, 4]

    # input g2
    tag2 = data[:, 9]

    # snr
    snr = data[:, 19]
    # snr = data[:, 20]
    idx0 = snr != 0
    # ssnr = snr[idx0]
    # idxs = ssnr >= snr_cut_s
    # idxe = ssnr <= snr_cut_e

    # the first 4 rows are the ellipticity,
    # the second 4 rows are the corresponding error bar,
    # the third 4 rows are the corresponding number of samples.
    res_arr1 = numpy.zeros((12, len(fg1)))
    res_arr2 = numpy.zeros((12, len(fg2)))

    for i in range(len(fg1)):
        idx11 = tag1 > fg1[i] - 0.0001
        idx12 = tag1 < fg1[i] + 0.0001
        ssnr = snr[idx11 & idx12 & idx0]
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
                G1 = data[:, 3]
                N1 = data[:, 10]
                U1 = data[:, 11]
                # V1 = data[:, 13]
                # V1.shape = (len(V1), 1)

                G1 = G1[idx11&idx12&idx0]
                N1 = N1[idx11&idx12&idx0]
                U1 = U1[idx11&idx12&idx0]
                G11 = numpy.append(G1[idxs], G1[idxe])
                N11 = numpy.append(N1[idxs], N1[idxe])
                U11 = numpy.append(U1[idxs], U1[idxe])
                weight1 = 1#snr[idx11&idx12&idxs&idxe]
                num1 = len(G11)
                # g1_h, g1_h_sig = Fourier_Quad().fmin_g(G1, N1, U1, mode=1, bin_num=8, sample=500)
                sig1 = []
                for k in range(20):
                    choice1 = numpy.random.randint(0, num1, 200000)
                    sig1.append(numpy.sum(G11[choice1]) / numpy.sum(N11[choice1]))
                g1_h_sig = numpy.std(sig1)
                g1_h = numpy.sum(G11 * weight1)/numpy.sum(N11 * weight1)
                # g1_h_sig = numpy.std(G1 * weight1/(N1 * weight1) - g1_h)/numpy.sqrt(num1)
            res_arr1[na, i] = g1_h
            res_arr1[na+4, i] = g1_h_sig
            res_arr1[na+8, i] = num1

    for i in range(len(fg2)):
        idx21 = tag2 > fg2[i] - 0.0001
        idx22 = tag2 < fg2[i] + 0.0001
        ssnr = snr[idx21 & idx22 & idx0]
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
                G2 = data[:, 8]
                N2 = data[:, 10]
                U2 = data[:, 11]
                # V1 = data[:, 13]
                # V1.shape = (len(V1), 1)
                G2 = G2[idx21&idx22&idx0]
                N2 = N2[idx21&idx22&idx0]
                U2 = U2[idx21&idx22&idx0]
                G22 = numpy.append(G2[idxs], G2[idxe])
                N22 = numpy.append(N2[idxs], N2[idxe])
                U22 = numpy.append(U2[idxs], U2[idxe])
                weight2 = 1#snr[idx21&idx22&idxs&idxe]
                num2 = len(G22)
                # g2_h, g2_h_sig = Fourier_Quad().fmin_g(G2, N2, U2, mode=2, bin_num=8, sample=500)
                g2_h = numpy.sum(G22 * weight2)/numpy.sum(N22 * weight2)
                sig2 = []
                for k in range(20):
                    choice2 = numpy.random.randint(0, num2, 200000)
                    sig2.append(numpy.sum(G22[choice2]) / numpy.sum(N22[choice2]))
                g2_h_sig = numpy.std(sig2)
                # g2_h_sig = numpy.std((G2 * weight2)/(N2 * weight2) - g2_h)/numpy.sqrt(num2)
            res_arr2[na, i] = g2_h
            res_arr2[na+4, i] = g2_h_sig
            res_arr2[na+8, i] = num2

    final_cache_path = path+'final_cache'
    numpy.savez(final_cache_path, res_arr1, res_arr2)

else:
    text = numpy.load(path+'final_cache.npz')
    res_arr1 = text['arr_0']
    res_arr2 = text['arr_1']
tm =time.time()

# fit the line
print('done\nbegin to plot the lines')

name = ['KSB', 'BJ02', 'Re-Gaussianization', 'Fourier_Quad']

mc_data_path = path + 'mc_data.xlsx'
mc_data = numpy.zeros((32, 1))
for i in range(3, 4):
    # Y = A*X ,   y = m*x+c
    # Y = [y1,y2,y3,...].T  the measured data
    # A = [[1,1,1,1,...]
    #         [y1,y2,y3..]].T
    # X = [c,m].T
    # C = diag[sig1^2, sig2^2, sig3^2, .....]
    # the inverse of C is used as weight of data
    # X = [A.T*C^(-1)*A]^(-1) * [A.T*C^(-1) *Y]
    A1 = numpy.column_stack((numpy.ones_like(fg1.T), fg1.T))
    Y1 = res_arr1[i].T
    C1 = numpy.diag(res_arr1[i+4]**2)

    A2 = numpy.column_stack((numpy.ones_like(fg2.T), fg2.T))
    Y2 = res_arr2[i].T
    C2 = numpy.diag(res_arr2[i+4]**2)

    L1 = numpy.linalg.inv(numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), A1))
    R1 = numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), Y1)
    L2 = numpy.linalg.inv(numpy.dot(numpy.dot(A2.T, numpy.linalg.inv(C2)), A2))
    R2 = numpy.dot(numpy.dot(A2.T, numpy.linalg.inv(C2)), Y2)

    sig_m1 = numpy.sqrt(L1[1, 1])
    sig_c1 = numpy.sqrt(L1[0, 0])
    sig_m2 = numpy.sqrt(L2[1, 1])
    sig_c2 = numpy.sqrt(L2[0, 0])
    e1mc = numpy.dot(L1, R1)
    e2mc = numpy.dot(L2, R2)
    # plot g1 line
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
    ax.plot(fg1, e1mc[1] * fg1 + e1mc[0], label=name[i], color='red')
    ax.plot([-0.1, 0.1], [-0.1, 0.1], label='y=x', color='blue')
    ax.scatter(fg1, res_arr1[i, :], c='black')
    for j in range(len(fg1)):
        ax.text(fg1[j], res_arr1[i, j], str(round(res_arr1[i + 8, j] / 1000, 1)) + "K", color="red")

    ax.text(0.1, 0.85, 'm=' + str(round(e1mc[1] - 1, 5))+'$\pm$'+str(round(sig_m1, 5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.8, 'c=' + str(round(e1mc[0], 5))+'$\pm$'+str(round(sig_c1, 5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.75, str(snr_cut_s)+"$\leq$"+"S/N, S/N"+"$\leq$" + str(snr_cut_e), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g1', fontsize=20)
    plt.ylabel('Est  g1', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.07, 0.07)
    plt.xlim(-0.07, 0.07)
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
    ax.plot(fg2, e2mc[1] * fg2 + e2mc[0], label=name[i], color='red')
    ax.plot([-0.1, 0.1], [-0.1, 0.1], label='y=x', color='blue')
    ax.scatter(fg2, res_arr2[i, :], c='black')
    for j in range(len(fg2)):
        ax.text(fg2[j], res_arr2[i, j], str(round(res_arr2[i + 8, j] / 1000, 1)) + "K", color="red")
    ax.text(0.1, 0.85, 'm=' + str(round(e2mc[1] - 1, 5))+'$\pm$'+str(round(sig_m2, 5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.8, 'c=' + str(round(e2mc[0], 5))+'$\pm$'+str(round(sig_c2, 5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.1, 0.75, str(snr_cut_s)+"$\leq$"+"S/N, S/N"+"$\leq$" + str(snr_cut_e), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g2', fontsize=20)
    plt.ylabel('Est  g2', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.07, 0.07)
    plt.xlim(-0.07, 0.07)
    nm = pic_path + name[i] + ".png"
    plt.savefig(nm)
    print('plotted g2')

    # m1, m2
    mc_data[i] = e1mc[1]
    mc_data[i+16] = e2mc[1]
    # delta_m1, delta_m2
    mc_data[i+4] = sig_m1
    mc_data[i+20] = sig_m2
    # c1, c2
    mc_data[i+8] = e1mc[0]
    mc_data[i+24] = e2mc[0]
    # delta_c1, delta_c2
    mc_data[i+12] = sig_c1
    mc_data[i+28] = sig_c2

if filter_type != 'none':
    mc_col = [filter_type]# + '_' + snr_s + '_' + snr_e ][0]
    if os.path.exists(mc_data_path):
        df = pandas.read_excel(mc_data_path)
        df[mc_col[0]] = mc_data
        df.to_excel(mc_data_path)
    else:
        dex = ['Km1', 'Bm1', 'Rm1', 'Fm1', 'Kdm1', 'Bdm1', 'Rdm1', 'Fdm1', 'Kc1', 'Bc1', 'Rc1', 'Fc1', 'Kdc1', 'Bdc1', 'Rdc1', 'Fdc1','Km2', 'Bm2', 'Rm2', 'Fm2', 'Kdm2', 'Bdm2', 'Rdm2', 'Fdm2', 'Kc2', 'Bc2', 'Rc2', 'Fc2', 'Kdc2', 'Bdc2', 'Rdc2', 'Fdc2']
        mc_df = pandas.DataFrame(data=mc_data, index=dex, columns=mc_col)
        mc_df.to_excel(mc_data_path)
te = time.time()

print("Complete")
print(tm-ts, te-tm)
