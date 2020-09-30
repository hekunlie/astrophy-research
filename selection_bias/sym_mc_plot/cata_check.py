import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
import h5py
from plot_tool import Image_Plot
import tool_box
import os
import matplotlib.pyplot as plt
from Fourier_Quad import Fourier_Quad


data_path = "D:/result_pic_single_field/anamoly_field"
fnms = os.listdir("D:/result_pic_single_field/anamoly_field")
for fnm in fnms:
    field_tag = 106
    # field_name = data_path + "/data_%d.hdf5"%field_tag
    field_name = data_path + "/%s"%fnm

    h5f = h5py.File(field_name, "r")
    data = h5f["/data"][()]

    col_shift = 0

    nstar = data[:,col_shift+4]
    imax = data[:,col_shift+5]
    jmax = data[:,col_shift+6]
    gf1 = data[:, col_shift+14]
    gf2 = data[:, col_shift+15]

    idx1 = nstar >= 12
    idx2 = imax < 48
    idx3 = jmax < 48
    idx4 = numpy.abs(gf1) <= 0.005
    idx5 = numpy.abs(gf2) <= 0.005
    idx_ = idx1 & idx2 & idx3 & idx4 & idx5

    src_num = idx_.sum()

    mg1_all = data[:,16][idx_]
    mg2_all = data[:,17][idx_]
    mn_all = data[:,18][idx_]
    mu_all = data[:,19][idx_]

    print(src_num, "Gal")

    gf1 = data[:, col_shift+14][idx_]
    gf2 = data[:, col_shift+15][idx_]


    # set up gf bins
    bin_num1 = 32
    bin_num2 = 32
    gf1_bin = numpy.linspace(-0.005, 0.005, bin_num1+1,dtype=numpy.float32)
    gf2_bin = numpy.linspace(-0.005, 0.005, bin_num2+1,dtype=numpy.float32)
    gf1_pts = (gf1_bin[1:] + gf1_bin[:-1])/2
    gf2_pts = (gf2_bin[1:] + gf2_bin[:-1])/2

    fq = Fourier_Quad(12,124)

    for i in range(bin_num1):
        idx_11 = gf1 >= gf1_bin[i]
        idx_12 = gf1 < gf1_bin[i+1]
        idx = idx_11 & idx_12
        num1 = idx.sum()

        gh1, gh1_sig = fq.find_shear(mg2_all[idx], mn_all[idx] - mu_all[idx], 8)[:2]

        print("%.5f, %.5f(%.5f) %d"%(gf1_pts[i],gh1, gh1_sig,num1))

        idx_11 = gf2 >= gf2_bin[i]
        idx_12 = gf2 < gf2_bin[i+1]
        idx = idx_11 & idx_12
        num2 = idx.sum()

        gh2, gh2_sig = fq.find_shear(mg2_all[idx], mn_all[idx] - mu_all[idx], 8)[:2]
        print("%.5f, %.5f(%.5f) %d"%(gf2_pts[i],gh2, gh2_sig,num2))
