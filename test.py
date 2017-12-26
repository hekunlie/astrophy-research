from Fourier_Quad import Fourier_Quad
import numpy
import os
from astropy.io import fits
import tool_box
import matplotlib.pyplot as plt

fq = Fourier_Quad(52, 124321)
p = fq.ran_pos(45, 9)
psf = fq.cre_psf(4, "Moffat")
ppsf = fq.pow_spec(psf)
shear = numpy.load("E:/selection_bias/parameters/shear.npz")["arr_0"]
data = numpy.load("F:/debug/80pts/data_cache.npz")["arr_0"]

fg1 = numpy.sort(shear[0])
fg2 = numpy.sort(shear[1])
FG1 = data[:, 3]
FG2 = data[:, 8]
FN = data[:, 10]
FU = data[:, 11]
FV = data[:, 12]
tag1 = data[:, 4]

# input g2
tag2 = data[:, 9]
a = 0
b = 200000
res_arr1 = numpy.zeros((3, len(fg1)))
res_arr2 = numpy.zeros((3, len(fg2)))

for i in range(len(fg1)):
    idx11 = tag1 > fg1[i] - 0.000001
    idx12 = tag1 < fg1[i] + 0.000001
    G1 = FG1[idx11 & idx12][a:b]
    N1 = FN[idx11 & idx12][a:b]
    U1 = FU[idx11 & idx12][a:b]
    num1 = len(G1)
    g1_h, g1_h_sig = fq.fmin_g(G1, N1, U1, mode=1, bin_num=8, signal=fg1[i])
    res_arr1[0, i] = g1_h
    res_arr1[1, i] = g1_h_sig
    res_arr1[2, i] = num1
    print(num1, fg1[i], g1_h, g1_h_sig)

# for i in range(len(fg2)):
#     idx21 = tag2 > fg2[i] - 0.000001
#     idx22 = tag2 < fg2[i] + 0.000001
#     G2 = FG2[idx21 & idx22][a:b]
#     N2 = FN[idx21 & idx22][a:b]
#     U2 = FU[idx21 & idx22][a:b]
#     num2 = len(G2)
#     g2_h, g2_h_sig = fq.fmin_g(G2, N2, U2, mode=2, bin_num=8, signal=fg2[i])
#     res_arr2[0, i] = g2_h
#     res_arr2[1, i] = g2_h_sig
#     res_arr2[2, i] = num2
#     print(num2, fg2[i], g2_h, g2_h_sig)
#
# e1mc = tool_box.data_fit(fg1, res_arr2[0], res_arr2[1])
# e2mc = tool_box.data_fit(fg2, res_arr2[0], res_arr2[1])
# tool_box.mcplot(fg2, res_arr2, fg2, res_arr2, e1mc, e2mc, '0', '10000')













# for j in range(10):
#     g2 = shear[j]
#     path = "F:/%d/"%j
#     files = os.listdir(path)
#     print(files, j, g2)
#     path = path+files[0]
#     img = fits.open(path)[0].data
#     pool = fq.divide_stamps(img)
#     # for i in range(a, b):
#     noise = fq.noise(0, 380.86)
#     mg1, mg2, mn, mu = fq.shear_est(pool[1], ppsf, noise, F=True)[:4]
#     print(mg1, mg2)
#     # g2_h[0,j], g2_h[1,j] = fq.fmin_g(mg2, mn, mu, mode=1, bin_num=8)
#     # print(j, g2, g2_h[0,j], g2_h[1,j])
# # mc = tool_box.data_fit(shear[1], g2_h[0], g2_h[1])
# #
# # tool_box.mcplot(shear[1], g2_h, shear[1], g2_h, mc, mc, '0', '1000')