import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
import h5py
import tool_box
import matplotlib.pyplot as plt
from scipy import signal


x,y = numpy.mgrid[-100:100,-100:100]
fxy = 10 + x + 12*y + 1*x**2 - 0.02*x*y - 2*y**2
plt.imshow(fxy)
plt.show()

my, mx = numpy.mgrid[-1:2, -1:2]
print(y.shape)
y_fit, x_fit = my.flatten(), mx.flatten()
xy = [0, 0]
kernel1 = tool_box.fit_2d_kernel(x_fit, y_fit, xy, 2)
print(kernel1)

fxy_rec = signal.convolve2d(fxy, kernel1,"same","wrap")
plt.imshow(fxy_rec[3:197,3:197]-fxy[3:197,3:197])
plt.colorbar()
plt.show()


my, mx = numpy.mgrid[-2:3, -2:3]
print(y.shape)
y_fit, x_fit = my.flatten(), mx.flatten()
xy = [0, 0]
kernel2 = tool_box.fit_2d_kernel(x_fit, y_fit, xy, 2)
print(kernel2)

fxy_rec = signal.convolve2d(fxy, kernel2,"same","wrap")
plt.imshow(fxy_rec[5:195,5:195]-fxy[5:195,5:195])
plt.colorbar()
plt.show()


h5f = h5py.File("./kernel.hdf5","w")
h5f["/kernel1"] = numpy.float32(kernel1)
h5f["/kernel2"] = numpy.float32(kernel2)
h5f.close()