# -*- coding: utf-8 -*-
import matplotlib
matplotlib.use('Agg')
from sys import path
path.append('/home/hklee/codes/')
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
            if ".xlsx" in i:
                paths.append(path + i)

        # collect the data from the files and put into 'data_list'
        # 'data' is the final array that contains all the data
        tc1 = time.time()
        data_list = tool_box.classify(paths, 10)[0]
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

    # correction
    KSB_cor = data[:, 0]**2 + data[:, 5]**2
    KSB_cor.shape = (len(KSB_cor), 1)
    BJ_cor = data[:, 1]**2 + data[:, 6]**2
    BJ_cor.shape = (len(BJ_cor), 1)
    RG_cor = data[:, 2]**2 + data[:, 7]**2
    RG_cor.shape = (len(RG_cor), 1)
    cor = [KSB_cor, BJ_cor, RG_cor]

    # snr
    snr = data[:, 13]
    snr.shape = (len(snr), 1)
    idxs = snr >= snr_cut_s
    idxe = snr <= snr_cut_e

    # input g1
    tag1 = data[:, 4]
    tag1.shape = (len(tag1), 1)

    # input g2
    tag2 = data[:, 9]
    tag2.shape = (len(tag2), 1)

    # the first 4 rows are the ellipticity,
    # the second 4 rows are the corresponding error bar,
    # the third 4 rows are the corresponding number of samples.
    res_arr1 = numpy.zeros((12, len(fg1)))
    res_arr2 = numpy.zeros((12, len(fg2)))

    for i in range(len(fg1)):
        idx11 = tag1 > fg1[i] - 0.001
        idx12 = tag1 < fg1[i] + 0.001
        for na in range(4):
            if na != 3:
                ellip1 = data[:, na]
                ellip1.shape = (len(ellip1), 1)
                idx13 = ellip1 != -10
                idx14 = ellip1 < 1
                idx15 = ellip1 > -1
                # the measured e1
                e1 = ellip1[idx11&idx12&idx13&idx14&idx15&idxs&idxe]
                num1 = len(e1)
                if na==0:
                    g1_h = numpy.mean(e1)
                    g1_h_sig = numpy.std(e1)/numpy.sqrt(num1)
                else:
                    #measured_esq = cor[na][idx11&idx12&idx13&idxs&idxe]*(-1)+2
                    g1_h = numpy.mean(e1)#/numpy.mean(measured_esq)
                    g1_h_sig = numpy.std(e1)/numpy.sqrt(num1)#/measured_esq - g1_h)
            else:
                G1 = data[:, 3]
                G1.shape = (len(G1), 1)
                N1 = data[:, 10]
                N1.shape = (len(N1), 1)
                U1 = data[:, 11]
                U1.shape = (len(U1), 1)
                # V1 = data[:, 13]
                # V1.shape = (len(V1), 1)
                G1 = G1[idx11&idx12&idxs&idxe]
                N1 = N1[idx11&idx12&idxs&idxe]
                U1 = U1[idx11&idx12&idxs&idxe]
                num1 = len(G1)
                g1_h, g1_h_sig = Fourier_Quad().fmin_g(G1, N1, U1, mode=1, bin_num=8, sample=500)
                # g1_h = numpy.sum(G1)/numpy.sum(N1)
                # g1_h_sig = numpy.std(G1/N1-g1_h)/numpy.sqrt(num1)
            res_arr1[na, i] = g1_h
            res_arr1[na+4, i] = g1_h_sig
            res_arr1[na+8, i] = num1

    for i in range(len(fg2)):
        idx21 = tag2 > fg2[i] - 0.001
        idx22 = tag2 < fg2[i] + 0.001
        for na in range(4):
            if na != 3:
                ellip2 = data[:, na+5]
                ellip2.shape = (len(ellip2), 1)
                idx23 = ellip2 != -10
                idx24 = ellip2 < 1
                idx25 = ellip2 > -1
                # the measured e1
                e2 = ellip2[idx21&idx22&idx23&idx24&idx25&idxs&idxe]
                num2 = len(e2)
                if na==0:
                    g2_h = numpy.mean(e2)
                    g2_h_sig = numpy.std(e2)/numpy.sqrt(num2)
                else:
                    measured_esq = cor[na][idx21&idx22&idx23&idxs&idxe]
                    g2_h = numpy.mean(e2)#/numpy.mean((2 - measured_esq ))
                    g2_h_sig = numpy.std(e2)/numpy.sqrt(num2)#/measured_esq - g2_h)
            else:
                G2 = data[:, 8]
                G2.shape = (len(G2), 1)
                N2 = data[:, 10]
                N2.shape = (len(N2), 1)
                U2 = data[:, 11]
                U2.shape = (len(U2), 1)
                # V1 = data[:, 13]
                # V1.shape = (len(V1), 1)
                G2 = G2[idx21&idx22&idxs&idxe]
                N2 = N2[idx21&idx22&idxs&idxe]
                U2 = U2[idx21&idx22&idxs&idxe]
                num2 = len(G2)
                g2_h, g2_h_sig = Fourier_Quad().fmin_g(G2, N2, U2, mode=2, bin_num=8, sample=500)
                # g2_h = numpy.sum(G2)/numpy.sum(N2)
                # g2_h_sig = numpy.std(G2/N2-g2_h)/numpy.sqrt(num2)
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

name = ['KSB', 'BJ', 'REGAUSS', 'Fourier_Quad']

mc_data_path = path + 'mc_data.xlsx'
mc_data = numpy.zeros((12, 2))
for i in range(4):
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
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

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
    ax.text(0.1, 0.75, str(snr_cut_s)+"$\leq$"+"S/N"+"$\leq$" + str(snr_cut_e), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g1', fontsize=20)
    plt.ylabel('Est  g1', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.07, 0.07)
    plt.xlim(-0.07, 0.07)
    nm1 = pic_path + name[i] + "_g1.png"
    plt.savefig(nm1)
    print('plotted g1')

    #plot g2 line
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

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
    ax.text(0.1, 0.75, str(snr_cut_s)+"$\leq$"+"S/N"+"$\leq$" + str(snr_cut_e), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g2', fontsize=20)
    plt.ylabel('Est  g2', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    plt.ylim(-0.07, 0.07)
    plt.xlim(-0.07, 0.07)
    nm2 = pic_path + name[i] + "_g2.png"
    plt.savefig(nm2)
    print('plotted g2')

    # m1, m2
    mc_data[i] = e1mc[1], e2mc[1]
    # delta_m1, delta_m2
    mc_data[i+4] = sig_m1, sig_m2
    # c1, c2
    mc_data[i+8] = e1mc[0], e2mc[0]
    # delta_c1, delta_c2
    mc_data[i+12] = sig_c1, sig_c2

if filter_type is not 'none':
    mc_col1 = filter_type + '/' + snr_s + '~' + snr_e + '/g1'
    mc_col2 = filter_type + '/' + snr_s + '~' + snr_e + '/g2'
    if os.path.exists(mc_data_path):
        df = pandas.read_excel(mc_data_path)
        df[mc_col1] = mc_data[:, 0]
        df[mc_col2] = mc_data[:, 1]
        df.to_excel(mc_data_path)
    else:
        col = [mc_col1, mc_col2]
        dex = ['Km', 'Bm', 'Rm', 'Fm', 'Kdm', 'Bdm', 'Rdm', 'Fdm', 'Kc', 'Bc', 'Rc', 'Fc', 'Kdc', 'Bdc', 'Rdc', 'Fdc']
        mc_df = pandas.DataFrame(data=mc_data, index=dex, columns=col)
        mc_df.to_excel(mc_data_path)
te = time.time()

print("Complete")
print(tm-ts, te-tm)
