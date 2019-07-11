from sys import path
# path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import numpy
import h5py
from plot_tool import Image_Plot
import tool_box
import matplotlib.pyplot as plt
from numpy import fft



def NFW_profile(size, scale_radius, rho_s):
    cen = int(size/2)
    z, y, x = numpy.mgrid[-cen:cen, -cen:cen, -cen:cen]
    radius = numpy.sqrt(x**2+y**2+z**2)/scale_radius
    idx = radius == 0
    radius[idx] = 1
    density = rho_s/radius/(1 + radius)**2
    density[idx] = 0
    density[idx] = density[cen-1:cen+2][cen-1:cen+2][cen-1:cen+2].sum()/26
    return density


size = 40
cen = int(size/2)
sig = 4
z,y,x = numpy.mgrid[-cen:cen, -cen:cen, -cen:cen]

gauss = tool_box.gauss_profile(size, sig, cen, cen)
gauss_f = tool_box.image_fft(gauss)
gauss_f[cen, cen] = gauss_f[cen, cen]*1000
gauss_r = tool_box.image_ifft(gauss_f).real
img = Image_Plot()
img.subplots(2, 2)
img.set_style()
# print(i,y[i])
fig = img.axs[0][0].imshow(gauss)
img.figure.colorbar(fig, ax=img.axs[0][0])
fig = img.axs[0][1].imshow(numpy.abs(gauss_f))
img.figure.colorbar(fig, ax=img.axs[0][1])
fig = img.axs[1][0].imshow(gauss_r)
img.figure.colorbar(fig, ax=img.axs[1][0])
fig = img.axs[1][1].imshow(gauss_r - gauss)
img.figure.colorbar(fig, ax=img.axs[1][1])
img.show_img()
img.close_img()

fr = numpy.exp(-(x**2+y**2+z**2)/2/25)*100
fr_f = numpy.abs(numpy.fft.fftshift(fft.fftn(fr)))
print(fr_f.shape)

for i in range(cen-10, cen+10):
    img = Image_Plot()
    img.subplots(1, 2)
    img.set_style()
    # print(i,y[i])
    fig = img.axs[0][0].imshow(fr_f[i])
    img.figure.colorbar(fig, ax=img.axs[0][0])
    fig = img.axs[0][1].imshow(fr[i])
    img.figure.colorbar(fig, ax=img.axs[0][1])
    img.show_img()
    img.close_img()
