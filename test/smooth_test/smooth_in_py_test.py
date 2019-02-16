import os
from sys import path
path.append("E:/Github/astrophy-research/my_lib")
import time
import numpy
import tool_box
from Fourier_Quad import Fourier_Quad
import matplotlib.pyplot as plt


size = 50
my, mx = numpy.mgrid[0:size, 0:size]
fq = Fourier_Quad(size, 123)
gal = numpy.exp(-((mx-size/2)**2+(my-size/2)**2)/2/4)*1000 + numpy.random.normal(0, 20, size*size).reshape((size, size))

gal_p = fq.pow_spec(gal)
gal_log = numpy.log10(gal_p)
gal_log_fit = tool_box.smooth(gal_log, size)
gal_fit = 10**gal_log_fit


row, col = 32, 32

new = numpy.zeros_like(gal_log)
invs = numpy.zeros((25, 6))
for row in range(2, size-2):
    for col in range(2, size-2):
        tag = 0
        x = []
        y = []
        z = []
        pk = 0
        for i in range(-2, 3):
            for j in range(-2, 3):
                if tag not in [0, 4, 20, 24]:
                    if (i+row)*size + j+col != (size+1)*size/2:
                        x.append(j)
                        y.append(i)
                        z.append(gal_log[i+row, j+col])
                    else:
                        pk = tag
                tag += 1

        x = numpy.array(x)
        y = numpy.array(y)
        z = numpy.array(z)

        f = z.sum()
        fx = numpy.sum(z*x)
        fy = numpy.sum(z*y)
        fx2 = numpy.sum(z*x*x)
        fxy = numpy.sum(z*x*y)
        fy2 = numpy.sum(z*y*y)

        n = len(x)
        x1 = x.sum()
        x2 = numpy.sum(x*x)
        x3 = numpy.sum(x*x*x)
        x4 = numpy.sum(x*x*x*x)

        y1 = y.sum()
        y2 = numpy.sum(y*y)
        y3 = numpy.sum(y*y*y)
        y4 = numpy.sum(y*y*y*y)

        xy = numpy.sum(x*y)
        x2y = numpy.sum(x*x*y)
        xy2 = numpy.sum(x*y*y)
        x3y = numpy.sum(x*x*x*y)
        xy3 = numpy.sum(x*y*y*y)
        x2y2 = numpy.sum(x*x*y*y)

        cov = numpy.array([[n, x1, y1, x2, xy, y2],
                          [x1, x2, xy, x3, x2y, xy2],
                          [y1, xy, y2, x2y, xy2, y3],
                          [x2, x3, x2y, x4, x3y, x2y2],
                          [xy, x2y, xy2, x3y, x2y2, xy3],
                          [y2, xy2, y3, x2y2, xy3, y4]])
        inv_cov = numpy.linalg.inv(cov)
        invs[pk] = inv_cov[0]

        f_z = numpy.array([f,fx,fy,fx2, fxy,fy2])
        para = numpy.sum(tool_box.inv_cov()[pk]*f_z)
        print(para)
        new[row, col] = para
        # print(para, gal_log[row, col])
        arr = numpy.zeros((size, size))
        arr[int(size/2), int(size/2)] = 2
        # for ii in range(len(x)):
        #     arr[row+y[ii], col+x[ii]] = 1
        #
        # plt.imshow(arr)
        # plt.show()

print(invs)
print("S")
for i in range(25):
    s = ""
    for j in range(6):
        if j != 0:
                s += ", %.10f"%invs[i,j]
        else:
            s += "[%.10f" % invs[i, j]
    s += "],"
    print(s)

plt.subplot(221)
plt.imshow(gal_log[2:size-2, 2:size-2])

plt.subplot(222)
plt.imshow(gal_log_fit[2:size-2, 2:size-2])
plt.subplot(223)
plt.imshow(new[2:size-2, 2:size-2])
plt.subplot(224)
plt.imshow(new[2:size-2, 2:size-2] - gal_log_fit[2:size-2, 2:size-2])
plt.colorbar()
plt.show()