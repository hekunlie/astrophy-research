# -*- coding: utf-8 -*-
from sys import path
path.append("E:/Github/astrophy-research/")
import matplotlib.pyplot as plt
from Fourier_Quad import *
import numpy
import pandas

mc_data = pandas.read_excel('E:/selection_bias/new_data/Gauss_4_7x7/mc_data.xlsx').values[:,7:]
legends = ['[20, maximum]','[30, maximum]','[40, maximum]','[50, maximum]','[60, maximum]']
name = ['KSB', 'BJ', 'RG', 'FQ']
for i in range(4):
    fig = plt.figure(figsize=(20,20))
    plt.subplot(221)
    plt.ylabel('m')
    plt.xlabel('snr')
    plt.errorbar([20, 30, 40, 50, 60], mc_data[i]-1, mc_data[i+4])
    plt.title(name[i]+'_g1_m')
    plt.subplot(222)
    plt.ylabel('c')
    plt.xlabel('snr')
    plt.errorbar([20, 30, 40, 50, 60], mc_data[i+8], mc_data[i+12])
    plt.title(name[i]+'_g1_c')
    plt.subplot(223)
    plt.ylabel('m')
    plt.xlabel('snr')
    plt.errorbar([20, 30, 40, 50, 60], mc_data[i+16]-1, mc_data[i+20])
    plt.title(name[i]+'_g2_m')
    plt.subplot(224)
    plt.ylabel('c')
    plt.xlabel('snr')
    plt.errorbar([20, 30, 40, 50, 60], mc_data[i+24], mc_data[i+28])
    plt.title(name[i]+'_g2_c')
    plt.show()


