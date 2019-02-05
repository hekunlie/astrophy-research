import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
import numpy
from astropy.io import fits
import time
import h5py
import tool_box
import matplotlib.pyplot as plt


def fun_xy(xi, yi, params, xypows):
    vals = 0
    for i in range(terms):
        vals += params[i]*xi**(xypows[i][0])*yi**(xypows[i][1])
    return vals


order = int(argv[1])
terms = int((order + 1) * (order + 2) / 2)
img_num = int(argv[2])
rng = numpy.random.RandomState(123)
num = 200
x, y = numpy.mgrid[0:num, 0:num]
xf, yf = x.reshape((int(num*num),)), y.reshape((int(num*num),))

ch_num = int(argv[3])
ch = rng.choice(numpy.arange(0,len(xf)), int(ch_num), replace=False)
# x, y = numpy.linspace(0, 20,int(num*num)), numpy.linspace(-10, 10, int(num*num))
# xf, yf = x, y
# powers of each terms
pows = [(i-j, j) for i in range(order+1) for j in range(i+1)]

# the coefficients
coeffs = [3000., 1., -1.]
fxy_ori = 0.
if terms - 2 > 1:
    for i in range(1,terms - 2):
        coeffs.append(numpy.random.uniform(-3, 3, 1)[0]/i/2)

# the f(x,y)
for i in range(terms):
    fxy_ori += coeffs[i]*x**(pows[i][0])*y**(pows[i][1])

cof = ""
for k in range(terms):
    cof += "%g, " % coeffs[k]
print("The true coeff: ", cof)

print(fun_xy(0,0, coeffs, pows), fun_xy(55,20, coeffs, pows))
dx = 0
scales = 1
yscale = 1
for i in range(img_num):

    noise = numpy.random.normal(0, 20, int(num*num)).reshape((num, num))

    fxy = fxy_ori + noise

    fxyf = fxy.reshape((int(num*num),))

    h5f = h5py.File("/home/hkli/work/cpp/test/data/%d.hdf5"%i,"w")
    # xf = numpy.arange(0, 10)
    # yf = xf + 0.5
    # fxyf = xf - 5
    paras, pows, cov, fv = tool_box.fit_2d((xf[ch]+dx)*scales, (yf[ch]+dx)*scales, fxyf[ch]*yscale, order)

    print("The esimated: ",fun_xy((dx+0)*scales, (dx+0)*scales, paras, pows)/yscale,
          fun_xy((dx+55)*scales,(dx+20)*scales, paras, pows)/yscale)
    cov = numpy.array(cov)
    fv = numpy.array(fv)

    h5f["/x"] = xf[ch]
    h5f['/y'] = yf[ch]
    h5f['/fxy'] = fxyf[ch]
    h5f["/cov"] = cov
    h5f["/vect"] = fv
    h5f.close()
    cov_inv = numpy.linalg.inv(cov)
    if i >= 0:

        print(xf.sum(), yf.sum())
        a = ""
        print("Right")
        for j in range(terms):
            a += "%g, " % fv[j][0]
        print(a)

        print("\n")
        print("COV")
        for j in range(terms):
            s = ""
            for k in range(terms):
                s += "%g, " % cov[j][k]
            print(s)

        print("\n")
        print("COV_INV")
        for j in range(terms):
            s = ""
            for k in range(terms):
                s += "%g, " % cov_inv[j][k]
            print(s)
        #print("\n")
    fit_fxy = 0
    for j in range(terms):
        fit_fxy += paras[j,0]*x**(pows[j][0])*y**(pows[j][1])

    residuals = fxy - fit_fxy
    s = ""
    for j in range(terms):
        s += "%g, "%paras[j,0]
    print(cof)
    print(s)
    print("\n")
    plt.figure(figsize=(12,12))
    plt.subplot(221)
    plt.imshow(fxy)
    plt.colorbar()
    plt.subplot(222)
    plt.imshow(fit_fxy)
    plt.colorbar()
    plt.subplot(223)
    plt.imshow(residuals)
    plt.colorbar()
    plt.subplot(224)
    plt.hist(residuals.reshape(int(num*num), ), 30)
    plt.savefig("/home/hkli/work/cpp/test/pic/%d.png"%i)
    plt.close()
