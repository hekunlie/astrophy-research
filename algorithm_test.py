from Fourier_Quad import *
import numpy
import pandas
import matplotlib.pyplot as plt
from scipy import optimize
import time

res1 = numpy.zeros((40,3))
res2 = numpy.zeros((40,3))
for k in range(-20,20):
    while True:
        G2 = numpy.random.normal(loc=0,scale=0.4,size=10000000)
        if numpy.abs(numpy.mean(G2))<1.e-6:
            break
    print(numpy.mean(G2))
    G2=G2-10.*k/1000.
    t1 = time.clock()
    g_h,sig = Fourier_Quad().fmin_g(G2,10,0,1,8,method=2,sample=1000)
    t2=time.clock()
    g_h2, sig2 = Fourier_Quad().fmin_g(G2, 10, 0, 1, 8, method=1, sample=1000)
    t3 = time.clock()
    res1[k]= k/1000.,g_h,sig
    res2[k] = k / 1000., g_h2, sig2
    print(res1[k],t2-t1)
    print(res2[k],t3-t2)
#numpy.savez('algorithm_test.npz',res)