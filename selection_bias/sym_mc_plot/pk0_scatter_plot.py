import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
# path.append('%s/work/fourier_quad/' % my_home)
path.append('D:/GitHub/astrophy-research/mylib')
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
import numpy
import matplotlib
import matplotlib.ticker as mtick
import h5py
from astropy.io import fits


# the relation between Pk0 and magnitude

h5f = h5py.File("E:/mask_12.hdf5","r")
mask = h5f["/data"][()]
h5f.close()
idx = mask > 0

mask = mask[idx]
ch = numpy.random.choice(numpy.array([i for i in range(len(mask))]), 50000, replace=False).tolist()

h5f = h5py.File("E:/flux2_ex1_12.hdf5","r")
pk0 = h5f["/data"][()][idx][ch]/48/60
h5f.close()

h5f = h5py.File("E:/flux2_ex2_12.hdf5","r")
pk0_fit = h5f["/data"][()][idx][ch]/48/60
h5f.close()

h5f = h5py.File("E:/mag_true_12.hdf5","r")
mag_true = -h5f["/data"][()][idx][ch]
h5f.close()

h5f = h5py.File("E:/snr_sex_12.hdf5","r")
snr_sex = h5f["/data"][()][idx][ch]
h5f.close()

print(mag_true.min(), mag_true.max())
mag_min = mag_true.min()
mag_max = mag_true.max()

num = 500

mag_bin = numpy.linspace(mag_min, mag_max, num+1)
mags = (mag_bin[1:] + mag_bin[:num])/2
print(mags.shape)
pks = numpy.zeros((4,num))

for i in range(num):
    idx1 = mag_true >= mag_bin[i]
    idx2 = mag_true < mag_bin[i+1]


    data = snr_sex[idx1&idx2]
    pks[0,i], pks[1,i] = data.min(), data.max()

    data = snr_sex[idx1&idx2]
    pks[2,i], pks[3,i] = data.min(), data.max()


matplotlib.rcParams["font.family"] = "serif"

img = Image_Plot(legend_size=14, fig_x=6, fig_y=4,plt_line_width=2,axis_linewidth=2,xpad=0.2)
img.subplots(1,1)
img.axis_type(0,"major",tick_len=6, tick_width=1.2)
img.axis_type(0,"minor",tick_len=3, tick_width=1.2)
img.axis_type(1,"major",tick_len=6, tick_width=1.2)


labels = ["P$_{k0}$", "P$_{k0,fit}$", "MAX(P$_{k0}$,P$_{k0,fit}$)"]

# img.axs[0][0].scatter(mag[ch], pk0[idx][ch], s=5, color="C0",label=labels[0])
# img.axs[0][0].scatter(mag[ch], pk0_fit[idx][ch], s=5, color="C1",label=labels[1])
img.axs[0][0].fill_between(mags, pks[0], pks[1], alpha=0.8,color="C0",label=labels[0])
# img.axs[0][0].fill_between(mags, pks[2], pks[3], alpha=0.8,color="C1",label=labels[1])

ys = img.axs[0][0].set_ylim()
xs = img.axs[0][0].set_xlim()

img.axs[0][0].set_ylim(pks.min(), pks.max())
# img.axs[0][0].set_xlim(22.5, 25.02)
# img.axs[0][0].set_ylim(22, xs[1])
img.axs[0][0].set_yscale("log")
img.axs[0][0].set_xlabel("Magnitude",fontsize=img.xy_lb_size)
img.axs[0][0].set_ylabel("$P_{k0}$",fontsize=img.xy_lb_size)
img.axs[0][0].legend(loc="lower left", frameon=False, fontsize=img.legend_size)
# img.save_img("E:/pk_scatter.pdf")
img.show_img()


exit()

h5f = h5py.File("E:/mask_12.hdf5","r")
mask = h5f["/data"][()]
h5f.close()
idx = mask > 0

source_num = mask.shape[0]
mask = mask[idx]
# ch = numpy.random.choice(numpy.array([i for i in range(len(mask))]), 50000, replace=False).tolist()
gal_labels = numpy.arange(0,source_num)[idx]

h5f = h5py.File("E:/flux2_ex1_12.hdf5","r")
pk0 = h5f["/data"][()][idx]/48/60
h5f.close()

h5f = h5py.File("E:/flux2_ex2_12.hdf5","r")
pk0_fit = h5f["/data"][()][idx]/48/60
h5f.close()

h5f = h5py.File("E:/mag_true_12.hdf5","r")
mag_true = -h5f["/data"][()][idx]
h5f.close()

h5f = h5py.File("E:/mag_auto_12.hdf5","r")
mag_auto = h5f["/data"][()][idx]
h5f.close()

h5f = h5py.File("E:/snr_sex_12.hdf5","r")
sex_snr = h5f["/data"][()][idx]
h5f.close()


idx3 = pk0 < 2

idx1 = sex_snr > 30
idx2 = sex_snr < 35

idx = idx1&idx2&idx3

print(gal_labels[idx].shape)
print(gal_labels[idx][:20])
print(mag_auto[idx][:20])
print(mag_true[idx][:20])
print(sex_snr[idx][:20])
print(pk0[idx][:20])
img = Image_Plot()
img.subplots(1,1)
img.axs[0][0].hist(pk0[idx], 100)
img.show_img()

fits_img = fits.open("E:/gal_chip_0000.fits")[0].data
size = int(fits_img.shape[0]/100)
fq = Fourier_Quad(size,12)

gals = fq.segment(fits_img)[5676]

img = Image_Plot()
img.subplots(1,1)
fig = img.axs[0][0].imshow(gals)
img.figure.colorbar(fig,ax=img.axs[0][0])
img.show_img()