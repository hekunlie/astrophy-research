# import matplotlib
# matplotlib.use("Agg")
import numpy
import os
# my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
# path.append('%s/work/fourier_quad/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import time
from Fourier_Quad import Fourier_Quad
# # import galsim
import matplotlib.pyplot as plt
from astropy.io import fits
import tool_box
# from mpi4py import MPI
import h5py
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm
from numpy import fft


f = h5py.File("E:/total.hdf5")
data = f["/mc1"].value
print(data)
print(list(f.keys()))


exit(0)
def pow_spec(image):
    image_ps = fft.fft2(image)
    return image_ps

def inv_pow(image):
    image_ps = fft.ifft2(image)
    return image_ps

def err(message):
    err_message = "%s is wrong"%message
    raise ValueError(err_message)

def mag(snr, r):
    return -numpy.log10(snr*numpy.sqrt(numpy.pi*(r/0.187)**2)*120)/0.4 + 34.5358
r = numpy.linspace(0.2, 1.2, 100)
plt.plot(r, mag(5, r))
plt.show()
exit()
x, y = numpy.mgrid[24-15:24+15,24-15:24+15]

fq = Fourier_Quad(48, 123)
gals = fq.segment(img)
# gal = gals[115]
my, mx = numpy.mgrid[0:48,0:48]
gal = numpy.exp(-((my-24)**2+(mx-24)**2)/2./25.)*10#+\
# hdu = fits.PrimaryHDU(gal)
# hdu.writeto("E:/gal.fits",overwrite=True)
     #numpy.exp(-((my-20)**2+(mx-24)**2)/2/100)*100
# points = fq.ran_pos(num=100, radius=20, g=(0.01, 0))[1]
# gal = fq.convolve_psf(pos=points, psf_scale=4, flux=1000, psf="Moffat")# + noise


gal_p = numpy.log10(fq.pow_spec(gal))
# hdu = numpy.log10(gal_p)
# smooth_hdu = tool_box.smooth(hdu,48)
# # w = fits.PrimaryHDU(hdu)
# # w.writeto("E:/test.fits",overwrite=True)
# test = fits.open("E:/smooth.fits")[0].data
# plt.subplot(131)
# plt.imshow((test-smooth_hdu)/test)
# plt.colorbar()
# plt.subplot(132)
# plt.imshow(smooth_hdu)
# plt.colorbar()
# plt.subplot(133)
# plt.imshow(test)
# plt.colorbar()
# plt.show()
# exit()
gal_ps = tool_box.smooth(gal_p, 48)

peak_pow = numpy.max(gal_p)
fit_area = numpy.ones((5,5))*peak_pow
fit_x, fit_y = numpy.mgrid[22:27, 22:27]
rim = fq.border(1)
num = numpy.sum(rim)
fsnr = numpy.sqrt(10**gal_p[24,24]*num/numpy.sum(10**gal_p*rim))
fsnr_f = numpy.sqrt(10**gal_ps[24,24]*num/numpy.sum(10**gal_p*rim))
print(fit_area.shape, fit_x.shape, gal_p.shape, x.shape)
print(gal_ps.max()/gal_p.max(), fsnr, fsnr_f)

fig = plt.figure(figsize=(14,14))
ax1 = fig.add_subplot(231)
ax1.imshow(gal)
ax2 = fig.add_subplot(232)
ax2.set_title("Powerspectrum, SNR=%.2f"%fsnr)
ax2.imshow((gal_p - gal_ps)/gal_p)
ax3 = fig.add_subplot(233)
ax3.imshow(gal_ps)
ax3.set_title("Smoothed powerspectrum, SNR$_{fit}$=%.2f"%fsnr_f)
ax4 = fig.add_subplot(234)
ax4.imshow(gal_p)
# surf = ax4.plot_surface(x, y, gal_p[24-15:24+15,24-15:24+15], rstride=1, cstride=1, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# fig.colorbar(surf, shrink=0.5, aspect=10)
# ax4.plot_surface(fit_x, fit_y, fit_area)

ax5 = fig.add_subplot(235, projection='3d')
surf = ax5.plot_surface(x, y, gal_p[24-15:24+15,24-15:24+15], rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)
ax5.plot_surface(fit_x, fit_y, fit_area)
ax5.set_title("Peak: %g"%peak_pow)
ax6 = fig.add_subplot(236, projection="3d")
ax6.plot_surface(x, y, gal_ps[24-15:24+15,24-15:24+15], rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax6.set_title("Peak: %g"%gal_ps.max())
plt.suptitle("$\\frac{SNR}{SNR_{fit}}$ = %g"%(gal_p.max()/gal_ps.max()), fontsize=20)
plt.subplots_adjust(wspace=0.1,hspace=0.1)
plt.savefig("E:/2.png",bbox_inches='tight')
plt.show()


# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# cpus = comm.Get_size()
#
# data_path = "/mnt/ddnfs/data_users/hkli/simu_test_gal_bigger/result/data/data_2.0sig/data_%d.hdf5"%rank
# f = h5py.File(data_path)
# data = f["/data"].value
# f.close()
# snr_f = data[:, 2]
# snr = data[:, 3]
# mag = data[:, 6]
# idx1 = mag< 25
# idx2 = mag>24.5
# idx = data[:, 1] > 0
# plt.figure(figsize=(12, 6))
# plt.subplot(221)
# plt.scatter(mag[idx], snr[idx], s=0.3)
# plt.ylim(10**(-4),10**2.5)
# plt.yscale("log")
# plt.subplot(222)
# plt.scatter(mag[idx], snr_f[idx], s=0.3)
# plt.ylim(10**(-4),10**2.5)
# plt.yscale("log")
# plt.subplot(223)
# plt.scatter(snr[idx], snr_f[idx], s=0.3)
# # plt.hist(x=snr[idx],bins=40)
# plt.subplot(224)
# plt.hist(x=snr_f[idx], bins=40)
#
# plt.savefig("/home/hkli/work/test/snr_hist_%d.png"%rank)
# plt.close()
#




