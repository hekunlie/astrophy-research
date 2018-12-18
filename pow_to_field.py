import numpy
from sys import path
path.append("D:/Github/astrophy-research/my_lib/")
from Fourier_Quad import Fourier_Quad
import matplotlib.pyplot as plt
import tool_box
from numpy import fft

def pow_spec(image):
    image_ps = fft.fft2(image)
    return image_ps

def inv_pow(image):
    image_ps = fft.ifft2(image)
    return image_ps



size = 200
fq = Fourier_Quad(size, 123)
my, mx = numpy.mgrid[-size:size,0:size]
# print(my,mx)
pows = numpy.exp(-((my**2+(mx-size)**2)/2/10**2))# + numpy.exp(-(((my-200)**2+(mx-100)**2)/2/100**2))*2000
pow_inv = tool_box.mirror_arr(pows)
pows = numpy.column_stack((pows, pow_inv))

my, mx = numpy.mgrid[-size:size,-size:size]
pows = numpy.exp(-((my**2+mx**2)/2/3**2))*10000
plt.imshow(pows)
plt.show()
theta = numpy.random.random_sample(int(size*size*2)).reshape((size*2,size))*numpy.pi*2
theta_inv = tool_box.mirror_arr(theta)
theta = numpy.column_stack((theta, theta_inv))

pows = pows*numpy.cos(theta) + 1j*pows*numpy.sin(theta)
delta = fft.ifft2(fft.fftshift(pows))
plt.imshow(delta.real)
plt.show()
plt.hist(delta.real.flatten(), 100)
plt.show()
