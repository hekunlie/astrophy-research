# -*- coding: utf-8 -*-

import os
import numpy
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import tool_box
from Fourier_Quad import *
from sys import argv
import shelve


ts =time.clock()
snr= 'SNR>10'

fg1 = numpy.linspace(-0.005, 0.005, 11)
fg2 = numpy.linspace(-0.005, 0.005, 11)

#the lenght of the intervel
dfg1 = fg1[1]-fg1[0]
dfg2 = fg2[1]-fg2[0]

paths   = []

#where the result data file are placed
path  = '/lmc/selection_bias/result/data/'
#where the result figures will be created
pic_path = '/lmc/selection_bias/result/pic/'

# check the final result cache
exist = os.path.exists(path+'final_cache.npz')
if exist:
    # 0: to use the result cache data existed to plot the line and estimate the bias, 'm' and 'c'
    # 1: run the program to classify the data ( it could be skipped if this result cache exists) and estimate the shear
    print('0: use the result cache data existed\n1: overwrite it')
    comm = int(input('0/1?:'))
else :
    print('no result cache exists')

if not exist or comm==1:
    # check the classification cache
    dict_cache_exist = os.path.exists(path+'dict_cache')
    if dict_cache_exist:
        print('0: use the classification cache data existed\n1: overwrite it')
        dict_comm = int(input('0/1?:'))
    else:
        print('no classification cache')

    if not dict_cache_exist or dict_comm==1:
        print ('starting...')
        files = os.listdir(path)
        for i in files:
            if ".xlsx" in i:
                paths.append(path + i)

        g1 = {}
        g2 = {}
        # for FQ method
        fn1 = {}
        fn2 = {}
        fu1 ={}
        fv1={}
        fu2 ={}
        fv2={}
        # name: 0:KSB,1:BJ,2:RG,3:F_Q
        # dict g1 = {'0':..{fg1[i]:[...]},{fg1[i+1]:[...]}....,'1': ...{fg1[i]:[...]},{fg1[i+1]:[...]}...,'2':...,'3':....}
        for name in range(4):
            for i in fg1:
                g1.setdefault(name, {})[i] = []
                if name == 3:
                    fn1[i] = []
                    fu1[i] = []
                    fv1[i] = []
            for i in fg2:
                g2.setdefault(name, {})[i] = []
                if name==3:
                    fn2[i] = []
                    fu2[i] = []
                    fv2[i] = []
        # collect the data from the files and put into 'data_list'
        # elements in data-list are the data arrays
        data_list = tool_box.classify(paths,10)[0]
        for k in range(len(data_list)):
            # put the data into the corresponding list
            data = data_list[k]

            # field distortion fg1
            tag1 = data[:,5]
            tag1.shape = (len(tag1),1)

            # field distortion fg2
            tag2 = data[:,11]
            tag2.shape = (len(tag2),1)

            for i in range(len(fg1)):
                idx1 = tag1 == fg1[i]
                for na in range(4):
                    ellip1 = data[:,na]
                    ellip1.shape = (len(ellip1),1)
                    g1[na][fg1[i]].extend(numpy.ndarray.tolist(ellip1[idx1]))
                    if na ==3:
                        n1 = data[:,na+1]
                        n1.shape = (len(n1),1)
                        u1 = data[:,12]
                        u1.shape = (len(u1), 1)
                        v1 = data[:,13]
                        v1.shape = (len(v1), 1)
                        fn1[fg1[i]].extend(numpy.ndarray.tolist(n1[idx1]))
                        fu1[fg1[i]].extend(numpy.ndarray.tolist(u1[idx1]))
                        fv1[fg1[i]].extend(numpy.ndarray.tolist(v1[idx1]))

            for i in range(len(fg2)):
                idx2 = tag2 == fg2[i]
                for na in range(4):
                    ellip2 = data[:,na+6]
                    ellip2.shape = (len(ellip2),1)
                    g2[na][fg2[i]].extend(numpy.ndarray.tolist(ellip2[idx2]))
                    if na==3:
                        n2 = data[:,na+7]
                        n2.shape = (len(n2),1)
                        u2 = data[:,12]
                        u2.shape = (len(u2), 1)
                        v2 = data[:,13]
                        v2.shape = (len(v2), 1)
                        fn2[fg2[i]].extend(numpy.ndarray.tolist(n2[idx2]))
                        fu2[fg2[i]].extend(numpy.ndarray.tolist(u2[idx2]))
                        fv2[fg2[i]].extend(numpy.ndarray.tolist(v2[idx2]))
         # create the cache of the classification
        dict = [g1,g2,fn1,fn2,fu1,fu2,fv1,fv2]
        dict_name = ['g1','g2','fn1','fn2','fu1','fu2','fv1','fv2']
        dict_cache = shelve.open(path+'dict_cache')
        for i in range(len(dict_name)):
            dict_cache[dict_name[i]] = dict[i]
        dict_cache.close()
        print("Classification complete")


    else:
        #load the classification cache
        print('loading classification cache')
        dict_cache = shelve.open(path+'dict_cache')
        g1 = dict_cache['g1']
        g2 = dict_cache['g2']
        fn1 = dict_cache['fn1']
        fn2 = dict_cache['fn2']
        fu1 = dict_cache['fu1']
        fu2 = dict_cache['fu2']
        fv1 = dict_cache['fv1']
        fv2 = dict_cache['fv2']
        dict_cache.close()

    # the first 4 rows are the ellipticity,
    # the second 4 rows are the corresponding error bar,
    # the third 4 rows are the corresponding number of samples.
    res_arr1 = numpy.zeros((12, len(fg1)))
    res_arr2 = numpy.zeros((12, len(fg2)))

    print('calculating shears ')
    for i in range(4):
        for m in range(len(fg1)):
            if i != 3:
                #for KSB, BJ, REGAUSS
                num1 = len(g1[i][fg1[m]])
                arr1 = numpy.array(g1[i][fg1[m]])
                ava1 = numpy.sum(arr1) / num1 / 1.6  # g1
                err1 = numpy.std(arr1) / numpy.sqrt(num1)   #error bar
                res_arr1[i, m] = ava1
                res_arr1[i + 4, m] = err1
                res_arr1[i + 8, m] = num1
            else:
                #for Fourier_Quad
                G1 = numpy.array(g1[i][fg1[m]])
                N1 = numpy.array(fn1[fg1[m]])
                U1 = numpy.array(fu1[fg1[m]])
                num1 = len(G1)
                bin_num =6
                g1_h,xi1_sq_min = Fourier_Quad().fmin_g(G1,U1,N1,mode=1,bin_num=bin_num)
                g1_h_p = Fourier_Quad().fmin_g(G1,U1,N1,mode=1,bin_num=bin_num,twi=2*xi1_sq_min,left=g1_h)[0]
                res_arr1[i, m] = g1_h
                res_arr1[i + 4, m] = g1_h_p-g1_h
                res_arr1[i + 8, m] = num1

        for m in range(len(fg2)):
            if i != 3:
                num2 = len(g2[i][fg2[m]])
                arr2 = numpy.array(g2[i][fg2[m]])
                ava2 = numpy.sum(arr2) / num2 / 1.6  # g2
                err2 = numpy.std(arr2) / numpy.sqrt(num2 )  #error bar
                res_arr2[i, m] = ava2
                res_arr2[i + 4, m] = err2
                res_arr2[i + 8, m] = num2
            else:
                G2 = numpy.array(g2[i][fg2[m]])
                N2 = numpy.array(fn2[fg2[m]])
                U2 = numpy.array(fu2[fg2[m]])
                num2 = len(G2)
                bin_num = 6
                g2_h,xi2_sq_min = Fourier_Quad().fmin_g(G2, U2, N2, mode=2,bin_num=bin_num)
                g2_h_p  = Fourier_Quad().fmin_g(G2,U2,N2,mode=2,bin_num=bin_num,twi=2*xi2_sq_min,left=g2_h)[0]
                res_arr2[i, m] = g2_h
                res_arr2[i + 4, m] = g2_h_p-g2_h
                res_arr2[i + 8, m] = num2

    final_cache_path = path+'final_cache'
    numpy.savez(final_cache_path,res_arr1,res_arr2)

else:
    text = numpy.load(path+'final_cache.npz')
    arr1 = text['arr_0']
    arr2 = text['arr_1']
tm =time.clock()
# fit the line
start1=0
end1 =0
start2=0
end2=0
print('done\nbegin to plot the lines')

fgn1 = fg1
fgn2 = fg2

name = ['KSB', 'BJ', 'REGAUSS', 'Fourier_Quad']
for i in range(4):
    # Y = A*X ,   y = m*x+c
    # Y = [y1,y2,y3,...].T  the measured data
    # A = [[1,1,1,1,...]
    #         [y1,y2,y3..]].T
    # X = [c,m].T
    # C = diag[sig1^2, sig2^2, sig3^2, .....]
    # the inverse of C is used as weight of data
    # X = [A.T*C^(-1)*A]^(-1) * [A.T*C^(-1) *Y]
    A1 = numpy.column_stack((numpy.ones_like(fgn1.T),fgn1.T))
    Y1  = arr1[i].T
    C1  = numpy.diag((arr1[i+4].T**2))
    A2 = numpy.column_stack((numpy.ones_like(fgn2.T),fgn2.T))
    Y2  = arr2[i].T
    C2  = numpy.diag((arr2[i+4].T**2))

    L1 = numpy.linalg.inv(numpy.dot(numpy.dot(A1.T,numpy.linalg.inv(C1)),A1))
    R1 = numpy.dot(numpy.dot(A1.T,numpy.linalg.inv(C1)),Y1)
    L2 = numpy.linalg.inv(numpy.dot(numpy.dot(A2.T,numpy.linalg.inv(C2)),A2))
    R2 = numpy.dot(numpy.dot(A2.T,numpy.linalg.inv(C2)),Y2)

    sig_m1 = numpy.sqrt(L1[1,1])
    sig_c1   = numpy.sqrt(L1[0,0])
    sig_m2 = numpy.sqrt(L2[1,1])
    sig_c2   = numpy.sqrt(L2[0,0])
    e1mc = numpy.dot(L1,R1)
    e2mc = numpy.dot(L2,R2)
    # plot g1 line
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    xmajorLocator = MultipleLocator(0.005)
    xmajorFormatter = FormatStrFormatter('%1.3f')
    xminorLocator = MultipleLocator(0.001)

    ymajorLocator   = MultipleLocator(0.005)
    ymajorFormatter = FormatStrFormatter('%1.3f')
    yminorLocator   = MultipleLocator(0.001)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)

    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)

    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)

    ax.errorbar(fgn1, arr1[i, :], arr1[i + 4, :], ecolor='black', elinewidth='1', fmt='none',capsize=2)
    ax.plot(fgn1, e1mc[1] * fgn1 + e1mc[0], label=name[i], color='red')
    ax.plot(fgn1, fgn1, label='y=x', color='blue')
    ax.scatter(fg1, res_arr1[i, :], c='black')
    for j in range(len(fg1)):
        ax.text(fg1[j], res_arr1[i, j], str(round(res_arr1[i + 8, j] / 1000, 1)) + "K", color="red")

    ax.text(0.2, 0.9, 'm=' + str(round(e1mc[1] - 1, 5))+'$\pm$'+str(round(sig_m1,5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.2, 0.85, 'c=' + str(round(e1mc[0], 5))+'$\pm$'+str(round(sig_c1,5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.2, 0.8, snr, color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g1', fontsize=20)
    plt.ylabel('Est  g1', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    nm1 = pic_path + name[i] + "_g1.png"
    plt.savefig(nm1)
    print ('plotted g1')

    #plot g2 line
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    xmajorLocator = MultipleLocator(0.01)
    xmajorFormatter = FormatStrFormatter('%1.3f')
    xminorLocator = MultipleLocator(0.002)

    ymajorLocator   = MultipleLocator(0.01)
    ymajorFormatter = FormatStrFormatter('%1.3f')
    yminorLocator   = MultipleLocator(0.002)

    ax.xaxis.set_major_locator(xmajorLocator)
    ax.xaxis.set_major_formatter(xmajorFormatter)

    ax.yaxis.set_major_locator(ymajorLocator)
    ax.yaxis.set_major_formatter(ymajorFormatter)

    ax.xaxis.set_minor_locator(xminorLocator)
    ax.yaxis.set_minor_locator(yminorLocator)

    ax.errorbar(fgn2, arr2[i, :], arr2[i + 4, :], ecolor='black', elinewidth='1', fmt='none',capsize =2)
    ax.plot(fgn2, e2mc[1] * fgn2 + e2mc[0], label=name[i], color='red')
    ax.plot(fgn2, fgn2, label='y=x', color='blue')
    ax.scatter(fg2, res_arr2[i, :], c='black')
    for j in range(len(fg1)):
        ax.text(fg2[j], res_arr2[i, j], str(round(res_arr2[i + 8, j] / 1000, 1)) + "K", color="red")
    ax.text(0.2, 0.9, 'm=' + str(round(e2mc[1] - 1, 5))+'$\pm$'+str(round(sig_m2,5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.2, 0.85, 'c=' + str(round(e2mc[0], 5))+'$\pm$'+str(round(sig_c2,5)), color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.2, 0.8, snr, color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g2', fontsize=20)
    plt.ylabel('Est  g2', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    nm2 = pic_path + name[i] + "_g2.png"
    plt.savefig(nm2)
    print('plotted g2')
te = time.clock()
print ("Complete")
print(tm-ts,te-tm)
