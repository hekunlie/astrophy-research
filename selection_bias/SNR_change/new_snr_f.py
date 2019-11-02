import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
path.append("E:/Github/astrophy-research/my_lib")
from astropy.io import fits
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
import matplotlib.pyplot as plt


size1 = 32
size2 = 128
pts_n = 50
flux = tool_box.mag_to_flux(22)

fq_noise = Fourier_Quad(256, 1110)
noise = fq_noise.draw_noise(0, 60)
nsp = noise.shape

fq_pts1 = Fourier_Quad(size1, 1231)
fq_pts2 = Fourier_Quad(size2, 1231)
pts = fq_pts1.ran_pos(pts_n, 6)

gal1 = fq_pts1.convolve_psf(pts, 4, flux/pts_n)+noise[int(nsp[0]/2-size1/2):int(nsp[0]/2+size1/2),int(nsp[1]/2-size1/2):int(nsp[1]/2+size1/2)]
gal2 = fq_pts2.convolve_psf(pts, 4, flux/pts_n)+noise[int(nsp[0]/2-size2/2):int(nsp[0]/2+size2/2),int(nsp[1]/2-size2/2):int(nsp[1]/2+size2/2)]
print(gal1.shape, gal2.shape)
gal1_p = fq_pts1.pow_spec(gal1)
snr1 = fq_pts1.snr_f(gal1)
gal2_p = fq_pts2.pow_spec(gal2)
snr2 = fq_pts2.snr_f(gal2)

rim1 = fq_pts1.border(1)
rim2 = fq_pts2.border(1)
nlevel1 = numpy.sum(gal1_p*rim1)
nlevel2 = numpy.sum(gal2_p*rim2)
n1 = numpy.sum(rim1)
n2 = numpy.sum(rim2)

plt.imshow(gal1-gal2[int(64-size1/2):int(64+size1/2),int(64-size1/2):int(64+size1/2)])
plt.show()
plt.close()
print(snr1, snr2)
print("%g, %g"%(nlevel1/n1, nlevel2/n2))
print(n1, n2)
plt.subplot(221)
plt.imshow(gal1_p)
plt.subplot(222)
plt.imshow(gal2_p)
plt.subplot(223)
plt.imshow(gal1)
plt.subplot(224)
plt.imshow(gal2)
plt.show()