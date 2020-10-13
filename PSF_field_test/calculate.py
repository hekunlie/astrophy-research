from sys import path
path.append('D:/GitHub/astrophy-research/mylib')
import numpy
from plot_tool import Image_Plot
import h5py
from Fourier_Quad import Fourier_Quad


h5f = h5py.File("D:/data_noise_free_0.hdf5","r")
mg1_0 = h5f["/mg1"][()]
mg2_0 = h5f["/mg2"][()]
mn_0 = h5f["/mn"][()]
mu_0 = h5f["/mu"][()]
h5f.close()

h5f = h5py.File("D:/data_noise_free_psf_r_0.hdf5","r")
mg1_1 = h5f["/mg1"][()]
mg2_1 = h5f["/mg2"][()]
mn_1 = h5f["/mn"][()]
mu_1 = h5f["/mu"][()]
h5f.close()



fq = Fourier_Quad(12,124)

gh1, gh1_sig = fq.find_shear(mg1_0, mn_0+mu_0, 8)[:2]
gh2, gh2_sig = fq.find_shear(mg2_0, mn_0-mu_0, 8)[:2]
print(gh1, gh1_sig)
print(gh2, gh2_sig)

gh1, gh1_sig = fq.find_shear(mg1_1, mn_1+mu_1, 8)[:2]
gh2, gh2_sig = fq.find_shear(mg2_1, mn_1-mu_1, 8)[:2]
print(gh1, gh1_sig)
print(gh2, gh2_sig)