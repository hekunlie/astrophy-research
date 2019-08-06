import numpy
import matplotlib.pyplot as plt
from numpy import fft
import scipy
from scipy.optimize import least_squares
from sys import path
path.append("E:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import tool_box
import h5py

lensing_low = numpy.loadtxt("E:/Github/astrophy-research/galaxy-galaxy lensing/GGL_calculation/lensing_low/data.dat")

cata_path = "F:/works/gg_lensing/cmass/fourier/old_cata/Tomo/"
h5f = h5py.File(cata_path+ "w_1_1.hdf5","r")
redshift = h5f["/Z"].value
h5f.close()
print(redshift.min(), redshift.max())
fq = Fourier_Quad(12, 124)
filen = 9
result = numpy.zeros((6, filen))
for i in range(filen):
    h5f = h5py.File(cata_path + "radius_%d.hdf5"%i,"r")
    data = h5f["/pair_data_0"].value
    h5f.close()
    crit_integ = data[:,4]*554.682135528
    mgt = data[:,0]*crit_integ
    mgx = data[:,1]*crit_integ
    mnut = data[:,2]
    mnux = data[:,3]
    print(data[:,-1].min())
    print(mgt.shape[0])
    img = Image_Plot()
    img.subplots(1, 2)
    gt, gt_sig = fq.fmin_g_new(mgt, mnut, 8, left=-180, right=180, fig_ax=img.axs[0][0])[:2]
    gx, gx_sig = fq.fmin_g_new(mgx, mnux, 8, left=-180, right=180, fig_ax=img.axs[0][1])[:2]
    # print(result[:,i])
    img.show_img()
    result[0,i] = gt
    result[1,i] = gt_sig
    result[2,i] = gx
    result[3,i] = gx_sig
    result[4,i] = data[:,-2].mean()


img = Image_Plot()
img.subplots(1,1)

names = ["old cata", "new cata",  "lensing low"]

img.axs[0][0].errorbar(result[4], result[0], result[1],capsize=3, c="C1",label="T")
img.axs[0][0].errorbar(result[4], result[2], result[3],capsize=3, c="C2",label="x")
img.axs[0][0].errorbar(lensing_low[:,0], lensing_low[:,1], lensing_low[:,2], capsize=3,c="C4",label="%s"%names[2])
img.axs[0][0].legend()
img.axs[0][0].set_yscale("log")
img.axs[0][0].set_xscale("log")
img.axs[0][0].set_ylim(0.1,180)
img.show_img()