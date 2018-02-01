import numpy
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import time
from Fourier_Quad import Fourier_Quad
from sys import argv
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def ellip(e,radius, theta):
    b = numpy.sqrt((1-e)/(1+e))*radius
    num = 0
    xy = numpy.zeros((2,100))
    while True:
        if num == 100:
            break
        y = numpy.random.normal(0, b,1)
        x = numpy.random.normal(0, radius,1)
        if x*x/radius/radius + y*y/b/b <=1:
            xy[0,num] = x
            xy[1, num] = y
            num += 1
    rotate = numpy.array([[numpy.cos(theta), -numpy.sin(theta)], [numpy.sin(theta), numpy.cos(theta)]])
    return numpy.dot(rotate, xy)

my,mx = numpy.mgrid[0:100,0:100]
arr = numpy.exp(-(my-50)**2/200-(mx-50)**2/400)
fig = plt.figure()
ax1 = Axes3D(fig)

ax1.plot_surface(mx, my, arr)
plt.show()
for i in range(10):
    xy = ellip(0.3, 10, numpy.pi/2)
    plt.scatter(xy[0], xy[1])
    plt.xlim(-12,12)
    plt.ylim(-12,12)
    plt.show()
    plt.close()