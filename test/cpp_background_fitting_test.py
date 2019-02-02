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

order = int(argv[1])
terms = int((order + 1) * (order + 2) / 2)

num = 200
x, y = numpy.mgrid[-num/2:num/2, -num/2:num/2]/10
xf, yf = x.reshape((int(num*num),)), y.reshape((int(num*num),))
pows = [(i-j, j) for i in range(order+1) for j in range(i+1)]
coeffs = [3000., 1., -1.]
fxy_ori = 0.
if terms - 2 > 1:
    for i in range(1,terms - 2):
        coeffs.append(numpy.random.uniform(-5,5,1)[0]/i/2)
for i in range(terms):
    fxy_ori += coeffs[i]*x**(pows[i][0])*y**(pows[i][1])
cof = ""
for k in range(terms):
    cof += "%g, " % coeffs[k]
print("The true coeff: ", cof)
for i in range(int(argv[2])):

    noise = numpy.random.normal(0, 10, int(num*num)).reshape((num, num))

    fxy = fxy_ori + noise

    fxyf = fxy.reshape((int(num*num),))

    h5f = h5py.File("/home/hkli/work/cpp/test/data/%d.hdf5"%i,"w")
    h5f["/x"] = xf
    h5f['/y'] = yf
    h5f['/fxy'] = fxyf
    h5f.close()

    paras, pows, cov, fv = tool_box.fit_2d(xf, yf, fxyf, order)
    cov_inv = numpy.linalg.inv(numpy.array(cov))
    if i >= 0:

        print(xf.sum(), yf.sum())
        a = ""
        for j in range(terms):
            a += "%g, " % fv[j][0]
        print(a)

        #print("\n")
        for j in range(terms):
            s = ""
            for k in range(terms):
                s += "%g, " % cov[j][k]
            #print(s)

        #print("\n")
        for j in range(terms):
            s = ""
            for k in range(terms):
                s += "%g, " % cov_inv[j][k]
            #print(s)
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
    plt.subplot(222)
    plt.imshow(fit_fxy)
    plt.subplot(223)
    plt.imshow(residuals)
    plt.colorbar()
    plt.subplot(224)
    plt.hist(residuals.reshape(40000, ), 30)
    plt.savefig("/home/hkli/work/cpp/test/pic/%d.png"%i)
    plt.close()