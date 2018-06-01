import numpy
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import time
# from Fourier_Quad import Fourier_Quad
from sys import argv
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
from scipy.optimize import curve_fit
# import tool_box
# from mpi4py import MPI
#
#

# show the defects of sextractor
a = numpy.loadtxt("F:/860599p_34.cat")
b = numpy.loadtxt("F:/test.cat")
minus = numpy.where(a[:,0] < 0)
print(a[minus],b[minus])
arr = fits.open("F:/860599p_34_source.fits")[0].data
arr = arr[1274-10:1274+10,1321-10:1321+10]
arr_ori = arr.copy()
arr[11,7] = 100
plt.subplot(121)
plt.title("SNR=-1.9")
plt.imshow(arr_ori)
plt.subplot(122)
plt.title("SNR=6.41")
plt.imshow(arr)
plt.show()
hdu = fits.PrimaryHDU(arr)
hdu.writeto("E:/a.fits",overwrite=True)


a = numpy.loadtxt("E:/821588p_13.cat")
arr = fits.open("E:/821588p_13_source.fits")[0].data
ab = numpy.where(a[:,0] > 1000)
x,y = int(a[ab][0,3]), int(a[ab][0,4])
print(a[ab])
print(x,y)
arr = arr[y-20:y+20,x-20:x+20]

arr_ori = arr.copy()
arr[19,13] = 0
plt.subplot(121)
plt.title("SNR=$10^{30}$")
plt.imshow(arr_ori)
plt.subplot(122)
plt.title("SNR=10.6")
plt.imshow(arr)
plt.show()
hdu = fits.PrimaryHDU(arr)
hdu.writeto("E:/a_big.fits",overwrite=True)

