# -*- coding: utf-8 -*-

import os
import numpy
import time
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy import optimize
from Fourier_Quad import *
from sys import argv

area,g1num,g2num= argv[1:4]

ts =time.clock()
snr= 'SNR>10'
g1num = int(g1num)
g2num = int(g2num)
fg1 = numpy.linspace(-0.005, 0.005, g1num)
fg2 = numpy.linspace(-0.01, 0.01, g2num)
dfg1 = fg1[1]-fg1[0] #the lenght of the intervel
dfg2 = fg2[1]-fg2[0]
paths   = []
path  = "/home/hklee/result/"+area+'/'           #where the result data file are placed
pic_path = '/home/hklee/result/pic/'+area+'/'    #where the result figures will be created
exist = os.path.exists(path+'cache.dat')

if exist:
    print ('0: use the cache data existed')
    print ('1: overwrite it')
    comm = int(input('0/1?:'))
else :
    print('no cache exists')

if not exist or comm==1:
    files = os.listdir(path)
    print ('starting...')
    for i in files:
        if ".txt" in i:
            paths.append(path + i)

    g1 = {}
    g2 = {}
    fn1 = {}  # for FQ method
    fn2 = {}
    fu1 ={}
    fv1={}
    fu2 ={}
    fv2={}
    # name: 0:KSB,1:BJ,2:RG,3:F_Q
    for name in range(4):  # dict g1/g2 = {'0':..{fg1/2[i]:[]},{fg1/2[i+1]:[]}..,'1':..,'2"...}
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

    for k in paths:  # put the data into the corresponding list
        if os.path.getsize(k)<1000:
            continue
        data = numpy.loadtxt(k, skiprows=1)[:,1:]
        tag1 = data[:,5]            #field distortion fg1
        tag1.shape = (len(tag1),1)
        tag2 = data[:,11]           #field distortion fg2
        tag2.shape = (len(tag2),1)

        for i in range(g1num):
            idx1 = tag1 < fg1[i] + dfg1/2
            idx2 = tag1 > fg1[i] - dfg1/2
            for na in range(4):
                ellip1 = data[:,na]
                ellip1.shape = (len(ellip1),1)
                g1[na][fg1[i]].extend(numpy.ndarray.tolist(ellip1[idx1&idx2]))
                if na ==3:
                    n1 = data[:,na+1]
                    n1.shape = (len(n1),1)
                    u1 = data[:,12]
                    u1.shape = (len(u1), 1)
                    v1 = data[:,13]
                    v1.shape = (len(v1), 1)
                    fn1[fg1[i]].extend(numpy.ndarray.tolist(n1[idx1&idx2]))
                    fu1[fg1[i]].extend(numpy.ndarray.tolist(u1[idx1&idx2]))
                    fv1[fg1[i]].extend(numpy.ndarray.tolist(v1[idx1&idx2]))
        for i in range(g2num):
            idx1 = tag2 < fg2[i] + dfg2/2
            idx2 = tag2 > fg2[i] - dfg2/2
            for na in range(4):
                ellip2 = data[:,na+6]
                ellip2.shape = (len(ellip2),1)
                g2[na][fg2[i]].extend(numpy.ndarray.tolist(ellip2[idx1&idx2]))
                if na==3:
                    n2 = data[:,na+7]
                    n2.shape = (len(n2),1)
                    u2 = data[:,12]
                    u2.shape = (len(u2), 1)
                    v2 = data[:,13]
                    v2.shape = (len(v2), 1)
                    fn2[fg2[i]].extend(numpy.ndarray.tolist(n2[idx1&idx2]))
                    fu2[fg2[i]].extend(numpy.ndarray.tolist(u2[idx1&idx2]))
                    fv2[fg2[i]].extend(numpy.ndarray.tolist(v2[idx1&idx2]))

    res_arr1 = numpy.zeros((12, g1num))    # the first 4 rows are the ellipticity,
    res_arr2 = numpy.zeros((12, g2num))    # the second 4 rows are the correspongding error bar,
                                            # the third 4 rows are the correspongding number of samples.
    for i in range(3,4):
        for m in range(len(fg1)):
            if i != 3:     #for KSB, BJ, REGAUSS
                num1 = len(g1[i][fg1[m]])
                arr1 = numpy.array(g1[i][fg1[m]])
                ava1 = numpy.sum(arr1) / num1 / 1.6  # g1
                err1 = numpy.std(arr1) / numpy.sqrt(num1)   #error bar
                res_arr1[i, m] = ava1
                res_arr1[i + 4, m] = err1
                res_arr1[i + 8, m] = num1
            else: #for Fourier_Quad
                G1 = numpy.array(g1[i][fg1[m]])
                num1 = len(G1)
                bin_num = 20
                inverse = range(int(bin_num/2-1),-1,-1)
                N1 = numpy.array(fn1[fg1[m]])
                U1 = numpy.array(fu1[fg1[m]])
                B1 = N1 + U1
                def fun(g):
                    G1_h = G1-B1*g
                    bins,num = Fourier_Quad().set_bins(G1_h,bin_num)[0:2]
                    n1 = num[0:int(bin_num/2)]
                    n2 = num[int(bin_num/2):][inverse]
                    return  numpy.sum((n1-n2)**2/(n1+n2))

                g1_h = optimize.fmin_cg(fun,[0.],disp=False)
                G1_h = G1-B1*g1_h
                num_in_bins= Fourier_Quad().set_bins(G1_h,bin_num)[1]
                n1 = num_in_bins[0:int(bin_num/2)]
                n2 = num_in_bins[int(bin_num/2):][inverse]
                sigma1 = 1./numpy.sqrt(2*numpy.sum((n1-n2)**2/(n1+n2)))
                # g1_h = numpy.mean(G1)/numpy.mean(N1)
                # sigma1 = numpy.std(G1)/numpy.mean(N1)/numpy.sqrt(num1)
                res_arr1[i, m] = g1_h
                res_arr1[i + 4, m] = sigma1
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
                num2 = len(G2)
                N2 = numpy.array(fn2[fg2[m]])
                U2 = numpy.array(fu2[fg2[m]])
                B2 =  N2 - U2
                def fun(g):
                    G2_h = G2 - B2 * g
                    bins, num = Fourier_Quad().set_bins(G2_h, bin_num)[0:2]
                    n1 = num[0:int(bin_num / 2)]
                    n2 = num[int(bin_num / 2):][inverse]
                    return numpy.sum((n1 - n2) ** 2 / (n1 + n2))

                g2_h = optimize.fmin_cg(fun, [0.], disp=False)
                G2_h = G2-B2*g2_h
                num_in_bins,size= Fourier_Quad().set_bins(G2_h,bin_num)[1:3]
                n1 = num_in_bins[0:int(bin_num/2)]
                n2 = num_in_bins[int(bin_num/2):][inverse]
                sigma2 = 1./numpy.sqrt(2*numpy.sum((n1-n2)**2/(n1+n2)))
                # g2_h = numpy.mean(G2)/numpy.mean(N2)
                # sigma2 = numpy.std(G2)/numpy.mean(N2)/numpy.sqrt(num2)
                res_arr2[i, m] = g2_h
                res_arr2[i + 4, m] = sigma2
                res_arr2[i + 8, m] = num2
    print ("Classification complete")

    with open(path+'cache.dat','wb') as cache:
        numpy.savetxt(cache,numpy.column_stack((res_arr1,res_arr2)))
else:
    text = numpy.loadtxt(path+'cache.dat')
    res_arr1 = text[:,0:g1num]
    res_arr2 = text[:,g1num:]
tm =time.clock()
# fit the line
start1=0
end1 =0
start2=0
end2=0

arr1 = res_arr1[:, start1:g1num-end1]
arr2 = res_arr2[:, start2:g2num-end2]
fgn1 = fg1[start1:g1num-end1]
fgn2 = fg2[start2:g2num-end2]

name = ['KSB', 'BJ', 'REGAUSS', 'Fourier_Quad']
for i in range(3,4):
    A1 = numpy.column_stack((numpy.ones_like(fgn1.T),fgn1.T))
    Y1  = arr1[i].T
    C1  = numpy.diag((arr1[i+4].T**2))  #sigma^2

    A2 = numpy.column_stack((numpy.ones_like(fgn2.T),fgn2.T))
    Y2  = arr2[i].T
    C2  = numpy.diag((arr2[i+4].T**2)) #sigma^2

    L1 = numpy.linalg.inv(numpy.dot(numpy.dot(A1.T,numpy.linalg.inv(C1)),A1))
    R1 = numpy.dot(numpy.dot(A1.T,numpy.linalg.inv(C1)),Y1)

    L2 = numpy.linalg.inv(numpy.dot(numpy.dot(A2.T,numpy.linalg.inv(C2)),A2))
    R2 = numpy.dot(numpy.dot(A2.T,numpy.linalg.inv(C2)),Y2)

    e1mc = numpy.dot(L1,R1)
    e2mc = numpy.dot(L2,R2)

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
    for j in range(g1num):
        ax.text(fg1[j], res_arr1[i, j], str(round(res_arr1[i + 8, j] / 1000, 1)) + "K", color="red")
    plt.xlim(-0.006, 0.006)
#    plt.ylim(-0.015, 0.015)
    m = 'm=' + str(round(e1mc[1] - 1, 5))
    c = 'c=' + str(round(e1mc[0], 5))
    ax.text(0.2, 0.9, m, color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.2, 0.85, c, color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.2, 0.8, snr, color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    plt.xlabel('True  g1', fontsize=20)
    plt.ylabel('Est  g1', fontsize=20)
    plt.title(name[i], fontsize=20)
    plt.legend(fontsize=15)
    nm1 = pic_path + name[i] + "_g1.png"
    plt.savefig(nm1)
    print ('plotted g1')

    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111)

    xmajorLocator = MultipleLocator(0.01)
    xmajorFormatter = FormatStrFormatter('%1.3f')
    xminorLocator = MultipleLocator(0.002)

    ymajorLocator   = MultipleLocator(0.01)
    ymajorFormatter = FormatStrFormatter('%1.3f')
    yminorLocator   = MultipleLocator(0.005)

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
    for j in range(g2num):
        ax.text(fg2[j], res_arr2[i, j], str(round(res_arr2[i + 8, j] / 1000, 1)) + "K", color="red")
    plt.xlim(-0.012, 0.012)
#    plt.ylim(-0.03, 0.03)
    m = 'm=' + str(round(e2mc[1] - 1, 5))
    c = 'c=' + str(round(e2mc[0], 5))

    ax.text(0.2, 0.9, m, color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
    ax.text(0.2, 0.85, c, color='green', ha='left', va='center', transform=ax.transAxes, fontsize=20)
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



                   
        
        
