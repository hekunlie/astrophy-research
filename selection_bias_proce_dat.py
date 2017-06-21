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
import shelve
import pandas

snr_cut = argv[1]
snr_cut = int(snr_cut)
g1num = 11
g2num = 11
fg1 = numpy.linspace(-0.005, 0.005, g1num)
fg2 = numpy.linspace(-0.005, 0.005, g2num)
dfg1 = fg1[1]-fg1[0] #the lenght of the intervel
dfg2 = fg2[1]-fg2[0]
paths   = []

results_path = "/lmc/selection_bias/result/data/"

exist = os.path.exists(path+'cache.dat')
if exist:#check the final result cache
    print('0: use the result cache data existed\n1: overwrite it')
    comm = int(input('0/1?:'))
else :
    print('no result cache exists')

if not exist or comm==1:

    dict_cache_exist = os.path.exists(path+'dict_cache.bak')  #check the classification cache
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
        fn1 = {}  # for FQ method
        fn2 = {}
        fu1 = {}
        fu2 = {}
        snr = []
        # name: 0:KSB,1:BJ,2:RG,3:F_Q
        for name in range(4):  # dict g1/g2 = {'0':..{fg1/2[i]:[]},{fg1/2[i+1]:[...the measured data...]}..,'1':..,'2"...}
            for i in fg1:
                g1.setdefault(name, {})[i] = []
                if name == 3:
                    fn1[i] = []
                    fu1[i] = []
            for i in fg2:
                g2.setdefault(name, {})[i] = []
                if name==3:
                    fn2[i] = []
                    fu2[i] = []

        for k in paths:  # put the data into the corresponding list
            # if os.path.getsize(k)<1000:
            #     continue
            data = pandas.read_excel(k).values
            tag1 = data[:,5]            #input g1
            tag1.shape = (len(tag1),1)
            tag2 = data[:,11]           #input g2
            tag2.shape = (len(tag2),1)
            snr_ori  = data[:,-1]
            idx_snr  = snr_ori > snr_cut
            for i in range(g1num):
                idx = tag1 = fg1[i]
                for na in range(4):
                    ellip1 = data[:,na]
                    ellip1.shape = (len(ellip1),1)
                    g1[na][fg1[i]].extend(numpy.ndarray.tolist(ellip1[idx&idx_snr]))
                    if na ==3:
                        n1 = data[:,na+1]
                        n1.shape = (len(n1),1)
                        u1 = data[:,12]
                        u1.shape = (len(u1), 1)
                        fn1[fg1[i]].extend(numpy.ndarray.tolist(n1[idx&idx_snr]))
                        fu1[fg1[i]].extend(numpy.ndarray.tolist(u1[idx&idx_snr]))
            for i in range(g2num):
                idx = tag2 = fg2[i]
                for na in range(4):
                    ellip2 = data[:,na+6]
                    ellip2.shape = (len(ellip2),1)
                    g2[na][fg2[i]].extend(numpy.ndarray.tolist(ellip2[idx&idx_snr]))
                    if na==3:
                        n2 = data[:,na+7]
                        n2.shape = (len(n2),1)
                        u2 = data[:,12]
                        u2.shape = (len(u2), 1)
                        fn2[fg2[i]].extend(numpy.ndarray.tolist(n2[idx&idx_snr]))
                        fu2[fg2[i]].extend(numpy.ndarray.tolist(u2[idx&idx_snr]))
         # create the cache of the classification
        dict = [g1,g2,fn1,fn2,fu1,fu2]
        dict_name = ['g1','g2','fn1','fn2','fu1','fu2']
        dict_cache = shelve.open(path+'dict_cache')
        for i in range(len(dict_name)):
            dict_cache[dict_name[i]] = dict[i]
        dict_cache.close()
        print("Classification complete")

    else:  #load the classification cache
        print('loading classification cache')
        dict_cache = shelve.open(path+'dict_cache')
        g1 = dict_cache['g1']
        g2 = dict_cache['g2']
        fn1 = dict_cache['fn1']
        fn2 = dict_cache['fn2']
        fu1 = dict_cache['fu1']
        fu2 = dict_cache['fu2']

    res_arr1 = numpy.zeros((12, g1num))    # the first 4 rows are the ellipticity,
    res_arr2 = numpy.zeros((12, g2num))    # the second 4 rows are the correspongding error bar,
                                            # the third 4 rows are the correspongding number of samples.
    print('calculating shears ')
    for i in range(4):
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
                bin_num =6
                inverse = range(int(bin_num/2-1),-1,-1)
                N1 = numpy.array(fn1[fg1[m]])
                U1 = numpy.array(fu1[fg1[m]])
                B1 = N1 + U1
                def G1_fun(g,twi):
                    G1_h = G1-B1*g
                    num = Fourier_Quad().set_bins(G1_h, bin_num, model=0)[1]
                    n1 = num[0:int(bin_num/2)]
                    n2 = num[int(bin_num/2):][inverse]
                    return  abs(numpy.sum((n1-n2)**2/(n1+n2))*0.5-twi)
                a = -0.1
                b = 0.1
                for i in range(4):
                    point = numpy.linspace[a, b, 10]
                    for k in range(len(point)):
                        mini0 = G1_fun(point[0],0)
                        kk    = 0
                        g1_h  = point[0]
                        mini = G1_fun(point[0],0)
                        if mini < mini0:
                            mini0 = mini
                            kk = k
                            g1_h = point[kk]
                    a = point[kk-1]
                    b = point[kk+1]
                a = g1_h
                b = 0.1
                twi = 2*mini0
                for i in range(4):
                    point = numpy.linspace[a, b, 10]
                    for k in range(len(point)):
                        mini0 = G1_fun(point[0],twi)
                        kk    = 0
                        del_g1  = point[0]
                        mini = G1_fun(point[k],twi)
                        if mini < mini0:
                            mini0 = mini
                            kk = k
                            del_g1 = point[kk]
                    a = point[kk-1]
                    b = point[kk+1]
                sigma1 =del_g1-g1_h

                # def fun(g,*args):
                #     half1,G,B,=args
                #     G1_h = G - B * g
                #     num = Fourier_Quad().set_bins(G1_h,bin_num,model=0)[1]
                #     n1 = num[0:int(bin_num/2)]
                #     n2 = num[int(bin_num/2):][inverse]
                #     return  abs(numpy.sum((n1-n2)**2/(n1+n2))*0.5-half1[0])
                # args=([0],G1,B1,)
                # g1_h = optimize.fmin(fun,[0.01],args=args,disp=False)[0]
                # args = ([fun(g1_h,[0],G1,B1,)*2], G1, B1,)
                # g1_hsig = optimize.fmin(fun,[0.01],args=args,disp=False)
                # print(g1_h,g1_hsig,abs(g1_h-g1_hsig))
                #sigma1 = abs(g1_h-optimize.fmin(fun,[0],args=args,disp=False)[0])
                # g1_h = numpy.mean(G1)/numpy.mean(N1)
                sigma1 = numpy.std(G1)/numpy.mean(N1)/numpy.sqrt(num1)
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
                def G1_fun(g,twi):
                    G2_h = G2-B2*g
                    num = Fourier_Quad().set_bins(G2_h, bin_num, model=0)[1]
                    n1 = num[0:int(bin_num/2)]
                    n2 = num[int(bin_num/2):][inverse]
                    return  abs(numpy.sum((n1-n2)**2/(n1+n2))*0.5-twi)
                a = -0.1
                b = 0.1
                for i in range(4):
                    point = numpy.linspace[a, b, 10]
                    for k in range(len(point)):
                        mini0 = G2_fun(point[0],0)
                        kk    = 0
                        g2_h  = point[0]
                        mini = G2_fun(point[0],0)
                        if mini < mini0:
                            mini0 = mini
                            kk = k
                            g2_h = point[kk]
                    a = point[kk-1]
                    b = point[kk+1]
                a = g2_h
                b = 0.1
                twi = 2*mini0
                for i in range(4):
                    point = numpy.linspace[a, b, 10]
                    for k in range(len(point)):
                        mini0 = G2_fun(point[0],twi)
                        kk    = 0
                        del_g1  = point[0]
                        mini = G2_fun(point[0],twi)
                        if mini < mini0:
                            mini0 = mini
                            kk = k
                            del_g1 = point[kk]
                    a = point[kk-1]
                    b = point[kk+1]
                sigma2 =del_g2-g2_h
                # def fun(g,*args):
                #     half2,G,B = args
                #     G2_h = G - B * g
                #     num = Fourier_Quad().set_bins(G2_h, bin_num,model=0)[1]
                #     n1 = num[0:int(bin_num / 2)]
                #     n2 = num[int(bin_num / 2):][inverse]
                #     return abs(numpy.sum((n1 - n2) ** 2 / (n1 + n2))*0.5-half2[0])
                # args = ([0],G2,B2)
                # g2_h = optimize.fmin(fun, [0.],args=args,disp=False)[0]
                # args = ([fun(g2_h,[0],G2,B2)*2], G2, B2,)
                # sigma2 = numpy.abs(g2_h-optimize.fmin(fun,[0],args=args,disp=False)[0])
                # g2_h = numpy.mean(G2)/numpy.mean(N2)
                #sigma2 = numpy.std(G2)/numpy.mean(N2)/numpy.sqrt(num2)
                res_arr2[i, m] = g2_h
                res_arr2[i + 4, m] = sigma2
                res_arr2[i + 8, m] = num2

    dict = [res_arr1,res_arr2]
    dict_name = ['arr1','arr2']
    dict_path = path+'/snr_cut/dict_cache_snr>%d'%snr_cut
    dict_cache = shelve.open(dict_path)
    if not os.paht.isdir(dict_path):
        os.makedirs(dict_path)
    for i in range(len(dict_name)):
        dict_cache[dict_name[i]] = dict[i]
    dict_cache.close()

else:
    dict_path = path + '/snr_cut/dict_cache_snr>%d' % snr_cut
    data_arr = numpy.load(dict_path)
    res_arr1 = data_arr["arr_0"]
    res_arr2 = data_arr["arr_0"]
tm =time.clock()
# fit the line
start1=0
end1 =0
start2=10
end2=10
print('done\nbegin to plot the lines')
arr1 = res_arr1[:, start1:g1num-end1]
arr2 = res_arr2[:, start2:g2num-end2]
fgn1 = fg1[start1:g1num-end1]
fgn2 = fg2[start2:g2num-end2]