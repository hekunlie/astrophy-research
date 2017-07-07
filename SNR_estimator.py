from Fourier_Quad import Fourier_Quad
import tool_box
import numpy
import os
from astropy.io import fits
import matplotlib.pyplot as plt
import copy
import scipy


signal_threshold = 576.38
img = fits.open('E:/cosmic shear/image/lsst_chips/gal_chip_99.fits')[0].data
#signal_threshold = tool_box.get_threshold(img,400,0.7,30,0)
pool = Fourier_Quad().divide_stamps(img,60)
mask = copy.copy(pool[8])
maxi = numpy.max(mask)
idx = mask <signal_threshold*1.5
mask[idx]=0
plt.imshow(mask)
plt.show()
arr_s = copy.copy(pool[8])
arr = numpy.zeros_like(mask)
collect = [1]
maxi = 0
for i in range(40):
    for j in range(40):
        signal = []
        if mask[i,j]!=0:
            p = tool_box.detect(mask,i,j,signal,0,(i,j))
            if p is not None and arr_s[p[0],p[1]]>maxi and len(signal)>4:
                maxi = arr_s[p[0],p[1]]
                collect[0]=(signal,p,maxi)
                for m in signal:
                    arr[m[0],m[1]]=arr_s[m[0],m[1]]
                    print(arr_s[m[0],m[1]],m)
                plt.subplot(121)
                plt.imshow(arr_s)
                plt.subplot(122)
                plt.imshow(arr)
                plt.show()
print(collect)
s=0
if collect[0] is not 1:
    for m in collect[0][0]:
        if arr_s[m[0],m[1]]>=collect[0][2]/2:
            s+=1
sig=numpy.sqrt(s/numpy.pi)/1.17
my,mx = numpy.mgrid[0:3,0:3]
my1,mx1 = numpy.mgrid[0:5,0:5]
g = 1./numpy.sqrt(2*numpy.pi)/sig*numpy.exp(-((my-1)**2+(mx-1)**2)/2/sig**2)
g1 = 1./numpy.sqrt(2*numpy.pi)/sig*numpy.exp(-((my1-2)**2+(mx1-2)**2)/2/sig**2)
det1=scipy.signal.convolve2d(arr_s,g,'same')
det2 =scipy.signal.convolve2d(arr_s,g1,'same')
pad = Fourier_Quad().border(15,60)
num = numpy.sum(pad)
noise1=numpy.sqrt(numpy.sum(det1**2*pad)/num)
noise2=numpy.sqrt(numpy.sum(det2**2*pad)/num)
print(det1[29,30]/noise1,det2[29,30]/noise2)
p1 =numpy.where(det1==numpy.max(det1))
p2 = numpy.where(det2==numpy.max(det2))
print(p1,p2)
plt.subplot(121)
plt.imshow(det1)
plt.subplot(122)
plt.imshow(det2)
plt.show()