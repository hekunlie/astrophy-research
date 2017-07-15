from Fourier_Quad import *
import numpy
import pandas
import matplotlib.pyplot as plt
from scipy import optimize
import time

res = numpy.zeros((40,3))
for k in range(-20,20):
    while True:
        G2 = numpy.random.normal(loc=0,scale=0.5,size=10000000)
        if numpy.abs(numpy.mean(G2))<1.e-6:
            break
    print(numpy.mean(G2))
    G2=G2-10.*k/100.
    g_h,sig = Fourier_Quad().fmin_g(G2,10,0,1,10)
    res[k]= k/100.,g_h,sig
    print(res[k])
numpy.savez('algorithm_test.npz',res)