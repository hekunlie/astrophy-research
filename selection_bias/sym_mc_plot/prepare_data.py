import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
import matplotlib
matplotlib.use("Agg")
from plot_tool import Image_Plot



total_path = argv[1]
cata_nm = "/total_raw.hdf5"

h5f = h5py.File(total_path + cata_nm,"r")
data = h5f["/data"][()]
h5f.close()

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

data_sub = data[idx_]

sigma = data_sub[:, col_shift+3]
hlr = data_sub[:,col_shift+7]
hla = data_sub[:,col_shift+8]
snr = hlr/numpy.sqrt(hla)
snr_sort = numpy.sort(snr)
snr_cut = [snr_sort[int(i*0.1*src_num)] for i in range(10)]

flux2 = data_sub[:, col_shift+10]
flux2_alt = data_sub[:, col_shift+11]
flux2_alt_sort = numpy.sort(flux2_alt)
flux2_alt_cut = [flux2_alt_sort[int(i*0.1*src_num)] for i in range(10)]

gf1 = data_sub[:, col_shift+14]
gf2 = data_sub[:, col_shift+15]

data_fq = data_sub[:,col_shift+16:]


bin_num = 20
gf_bin = numpy.linspace(-0.005, 0.005, bin_num+1,dtype=numpy.float32)
gf_pts = (gf_bin[1:] + gf_bin[:-1])/2

h5f = h5py.File(total_path + "/cutoff.hdf5","w")
h5f["/gf_bin"] = gf_pts
h5f["/snr"] = snr_cut
h5f["/flux2_alt"] = flux2_alt_cut

print(gf_pts)
print(h5f["/snr"][()])
print(h5f["/flux2_alt"][()])

for i in range(bin_num):
    idx_11 = gf1 >= gf_bin[i]
    idx_12 = gf1 < gf_bin[i+1]
    idx = idx_11 & idx_12
    num1 = idx.sum()

    h5f["/%d/mg_1"%i] = data_fq[:,0][idx]
    h5f["/%d/mn_1"%i] = data_fq[:,2][idx]
    h5f["/%d/mu_1"%i] = data_fq[:,3][idx]
    h5f["/%d/snr_1"%i] = snr[idx]
    h5f["/%d/flux2_alt_1"%i] = flux2_alt[idx]

    idx_11 = gf2 >= gf_bin[i]
    idx_12 = gf2 < gf_bin[i+1]
    idx = idx_11 & idx_12
    num2 = idx.sum()

    h5f["/%d/mg_2"%i] = data_fq[:,1][idx]
    h5f["/%d/mn_2"%i] = data_fq[:,2][idx]
    h5f["/%d/mu_2"%i] = data_fq[:,3][idx]
    h5f["/%d/snr_2"%i] = snr[idx]
    h5f["/%d/flux2_alt_2"%i] = flux2_alt[idx]

    print(i, gf_bin[i], gf_bin[i + 1], num1, num2)

h5f.close()