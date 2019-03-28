import numpy
import os
from sys import path,argv
path.append("E:/Github/astrophy-research/my_lib")
from Fourier_Quad import Fourier_Quad
import tool_box
import time
from astropy.io import fits
import matplotlib.pyplot as plt
import h5py


rng = numpy.random.RandomState(123)
fq = Fourier_Quad(12,122)

bin_num = 20
gh_num = 150
gh = numpy.linspace(-0.07, 0.07, gh_num)


signals = [0.05, -0.05]
sigmas = [2, 2]
nums = [50*1.e4, 50*1.e4]

datas = [rng.normal(signals[i], sigmas[i], int(nums[i])).reshape((int(nums[i]),1)) for i in range(len(signals))]

for i in range(len(datas)):
    if i == 0:
        data = datas[i]
    else:
        data = numpy.row_stack((data, datas[i]))
    print(data.shape)
print(bin_num, data.shape)
bins = fq.set_bin(data, bin_num)
print("Bin length: ", bins.shape)
plt.scatter(bins, [0 for i in range(len(bins))])
plt.show()


# each single signal
for i in range(len(datas)):
    chisq = []
    for j in range(gh_num):
        chisq.append(fq.G_bin(datas[i], 1, gh[j], bins, 0))
    plt.scatter(gh, chisq)
    plt.show()
    plt.close()
    est_g, est_g_sig = fq.fmin_g_new(datas[i], 1, bin_num)
    print(signals[i], est_g, est_g_sig)
chisq = []
for i in range(gh_num):
    chisq.append(fq.G_bin(data, 1, gh[i], bins, 0))
plt.figure(figsize=(16,12))
plt.scatter(gh, chisq)
plt.show()
plt.close()

