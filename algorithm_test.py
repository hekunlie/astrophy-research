from Fourier_Quad import *
import numpy
import pandas
import matplotlib.pyplot as plt
from scipy import optimize
import time

# res = numpy.zeros((20,4))
#
# for k in range(-10,10):
#     while True:
#         G2 = numpy.random.normal(loc=0,scale=0.4,size=1000000)
#         if numpy.abs(numpy.mean(G2))<1.e-5:
#             break
#     print(numpy.mean(G2))
#     G2=G2+10.*k/1000.
#     t1 = time.clock()
#     g_h,sig = Fourier_Quad().fmin_g(G2,10,0,1,8,method=2,sample=1000)
#     t2=time.clock()
#
#     res[k]= k/1000.,g_h,abs(g_h-k/1000.),sig
#     print(res[k],t2-t1)
# numpy.savez('algorithm_test.npz',res)

res = numpy.load('algorithm_test.npz')['arr_0']
mean = numpy.mean(res[:,2])*100

plt.scatter(res[:,0],res[:,2]*100,label = r'$\| residual \|$',linewidths=2,color='r')
#plt.yticks(numpy.linspace(0,0.012,5), ('0','$0.003*10^{-3}$','$0.006*10^{-3}$','$0.009*10^{-3}$','$0.012*10^{-3}$'))
plt.yticks(numpy.linspace(0,0.003,4), ('0','$0.001*10^{-2}$','$0.002*10^{-2}$','$0.003*10^{-2}$'))
plt.xlim(-0.011,0.011)
plt.ylim(-0.0001,0.003)
plt.xlabel('true signal',fontsize =15)
plt.ylabel(r'$\|residual\|$',fontsize =15)
plt.plot([-0.011,0.011],[numpy.mean(res[:,2])*100,numpy.mean(res[:,2])*100],label = r"$mean \approx 0.006*10^{-3}$",c='b')
plt.legend(fontsize = 13)
plt.show()