import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box
from astropy.io import fits
import os
import matplotlib.pyplot as plt


h5f = h5py.File("E:/data_cross_term_espf_0.hdf5","r")
data_ct = h5f["/data"][()]
h5f.close()

h5f = h5py.File("E:/data_noise_residual_espf_0.hdf5","r")
data_nr = h5f["/data"][()]
h5f.close()

fq = Fourier_Quad(12,124)
gh,chisq = fq.get_chisq_range(data_ct[:,0], data_ct[:,2]+)