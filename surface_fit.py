import matplotlib.pyplot as plt
import numpy
from astropy.io import fits
import copy
import tool_box
import time

order = 2
turns = int((order+1)*(order+2)/2)
power = numpy.zeros((2,turns))
k=0
for i in range(order+1):
    for j in range(i+1):
        power[0,k] = i-j
        power[1,k] = j
        k += 1
# print("xy power: \n", power)
my, mx = numpy.mgrid[-2:3, -2:3]
x, y = mx.reshape((1, 25)), my.reshape((1, 25))

# x = numpy.delete(x, [0, 4, 12, 20, 24])
# y = numpy.delete(y, [0, 4, 12, 20, 24])

set_zeros = [0, 4, 20, 24]
for i in set_zeros:
    x[0,i] = 0
    y[0,i] = 0
# print(y)
powx = numpy.zeros((turns, turns))
powy = numpy.zeros((turns, turns))
for m in range(turns):
    for n in range(turns):
        powx[m, n] = power[0, m] + power[0, n]
        powy[m, n] = power[1, m] + power[1, n]
# print(powx)
# print(powy)
print("x: \n", x)
print("y: \n", y)
fxy = numpy.zeros((turns, 25))
for i in range(turns):
    fxy[i] = x**power[0, i]*y**power[1, i]
print("fxy: \n",fxy)
for i in range(turns):
    s = ""
    for j in range(25):
        s += "%.1f, "%fxy[i,j]
    print("{ %s },"%s)
print("\n")
pts_list = []  # for checking

coeff_inv = numpy.zeros((25,turns))
for i in range(25):
    nx = copy.copy(x)
    ny = copy.copy(y)
    if i in [0, 4, 20, 24]:
        pts = 21
    else:
        nx[0, i] = 0
        ny[0, i] = 0
        pts = 20
    xy_arr = numpy.zeros((turns,turns))
    for m in range(turns):
        for n in range(turns):
            xy_arr[m,n] = numpy.sum(nx**(powx[m,n])*ny**(powy[m,n]))
    xy_arr[0,0] = pts

    pts_list.append(pts)

    inv_arr = numpy.linalg.inv(xy_arr)
    coeff_inv[i] = inv_arr[0]
    # print(coeff_inv)
t2 = time.clock()
for i in range(25):
    s = ""
    for j in range(turns):
        s += "%.8f, "%coeff_inv[i,j]
    print("{ %s }," % s)

# x4 = numpy.sum(x ** 4)
# x3y = numpy.sum(x ** 3 * y)
# x2y2 = numpy.sum(x ** 2 * y ** 2)
# xy3 = numpy.sum(x * y ** 3)
# y4 = numpy.sum(y ** 4)
# x3 = numpy.sum(x ** 3)
# x2y = numpy.sum(x ** 2 * y)
# x2 = numpy.sum(x ** 2)
# xy2 = numpy.sum(x * y ** 2)
# xy = numpy.sum(x * y)
# y3 = numpy.sum(y ** 3)
# y2 = numpy.sum(y ** 2)
# x1 = numpy.sum(x)
# y1 = numpy.sum(y)
# n = 20
#
# arr = numpy.array([[x4, x3y, x2y2, x3, x2y, x2], [x3y, x2y2, xy3, x2y, xy2, xy], [x2y2, xy3, y4, xy2, y3, y2],
#                    [x3, x2y, xy2, x2, xy, x1], [x2y, xy2, y3, xy, y2, y1], [x2, xy, y2, x1, y1, n]])
# print(arr)

size = 20
cen = int((size*size + size)/2)
for i in range(int(size/2)-3, int(size/2)+3):
    for j in range(int(size/2)-3, int(size/2)+3):
        arr = numpy.zeros((size, size))
        arr[int(size/2),int(size/2)] = 0.5
        pos = []
        tag = 0
        pk = 0
        for m in range(-2,3):
            p = (i + m + size) % size
            for n in range(-2,3):
                q = (j + n + size) % size

                if tag not in [0,4,20,24]:#abs(m) != 2 or abs(n) != 2:
                    if p*size+q != cen:
                        pos.append((p,q))
                    else:
                        pk = tag
                tag+=1
        print(pk,len(pos))
        for xy in pos:
            arr[xy[0], xy[1]] = 1
        plt.imshow(arr)
        plt.grid(which='minor')
        plt.show()
