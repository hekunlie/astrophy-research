from sys import path
#path.append('E:/Github/astrophy-research')
from Fourier_Quad import Fourier_Quad
import numpy
import pandas
import matplotlib.pyplot as plt
from scipy import optimize
import time
import os


mean_method = numpy.zeros((11,6))
smpdf = numpy.zeros((11,6))
files = os.listdir('E:/Github/astrophy-research/200K_gals_test_nonoise/')
path = []
for file in files:
    if 'gal_test_' in file:
        path.append(file)
print(path)
new = []
old = []
for i in range(11):
    print(i)
    data_matrix = numpy.load(path[i])['arr_0']
    G1 = data_matrix[:,0]
    G2 = data_matrix[:,1]
    N = data_matrix[:,2]
    U = data_matrix[:,3]
    g1 = data_matrix[0,-2]
    g2 = data_matrix[0, -1]
    g1_est_mean = numpy.sum(G1)/numpy.sum(N)
    g1_est_mean_sig = numpy.std(G1/N)/numpy.sqrt(len(N))
    g2_est_mean = numpy.sum(G2) / numpy.sum(N)
    g2_est_mean_sig = numpy.std(G2 / N) / numpy.sqrt(len(N))
    mean_method[i] = g1_est_mean,g1_est_mean_sig,g2_est_mean,g2_est_mean_sig,g1,g2
    t1=time.clock()
    g1_est,sig1  = Fourier_Quad().fmin_g(G1,N,U,mode=1,bin_num=8,sample=20000)[0:2]
    t2=time.clock()
    g1_est11, sig11 = Fourier_Quad().fmin_g(G1, N, U, mode=1, bin_num=8, method = 1,sample=20000)[0:2]
    t3=time.clock()
    g2_est,sig2 = Fourier_Quad().fmin_g(G2, N, U, mode=2, bin_num=8,sample=20000)[0:2]
    g2_est22, sig22 = Fourier_Quad().fmin_g(G2, N, U, mode=2, bin_num=8, method =1,sample=20000)[0:2]


    #print(g1,g1_est_mean,g2,g2_est_mean)
    print(g1,abs(g1-g1_est),abs(g1-g1_est11))
    print(g1, g1_est, g1_est11,t2-t1)
    print(g2, abs(g2 - g2_est), abs(g2 - g2_est22))
    print(g2, g2_est, g2_est22,t3-t2)
    old.append(g2-g2_est)
    new.append(g2-g2_est22)
    smpdf[i] = g1_est,sig1,g2_est,sig2,g1,g2
    mean_method[i] = g1_est_mean,g1_est_mean_sig,g2_est_mean,g2_est_mean_sig,g1,g2

print(numpy.std(old),numpy.std(new))
# numpy.savez('cache.npz',mean_method,smpdf)
#
# result = numpy.load('cache.npz')
# mean_method = result['arr_0']
# smpdf = result['arr_1']
# og1= mean_method[:,-2]
# og2 = mean_method[:,-1]
# mA1 = numpy.column_stack((numpy.ones_like(og1),og1))
# mY1  = mean_method[:,0]
# mC1  = numpy.diag(mean_method[:,1]**2)
# mA2 = numpy.column_stack((numpy.ones_like(smpdf[:,-1]),smpdf[:,-1]))
# mY2  = mean_method[:,2]
# mC2  = numpy.diag(mean_method[:,3]**2)
#
# mL1 = numpy.linalg.inv(numpy.dot(numpy.dot(mA1.T,numpy.linalg.inv(mC1)),mA1))
# mR1 = numpy.dot(numpy.dot(mA1.T,numpy.linalg.inv(mC1)),mY1)
# mL2 = numpy.linalg.inv(numpy.dot(numpy.dot(mA2.T,numpy.linalg.inv(mC2)),mA2))
# mR2 = numpy.dot(numpy.dot(mA2.T,numpy.linalg.inv(mC2)),mY2)
#
# msig_m1 = numpy.sqrt(mL1[1, 1])
# msig_c1 = numpy.sqrt(mL1[0, 0])
# msig_m2 = numpy.sqrt(mL2[1, 1])
# msig_c2 = numpy.sqrt(mL2[0, 0])
# me1mc = numpy.dot(mL1, mR1)
# me2mc = numpy.dot(mL2, mR2)
#
# pA1 = numpy.column_stack((numpy.ones_like(og1),og1))
# pY1  = smpdf[:,0]
# pC1  = numpy.diag(smpdf[:,1]**2)
# pA2 = numpy.column_stack((numpy.ones_like(smpdf[:,-1]),smpdf[:,-1]))
# pY2  = smpdf[:,2]
# pC2  = numpy.diag(smpdf[:,3]**2)
#
# pL1 = numpy.linalg.inv(numpy.dot(numpy.dot(pA1.T,numpy.linalg.inv(pC1)),pA1))
# pR1 = numpy.dot(numpy.dot(pA1.T,numpy.linalg.inv(pC1)),pY1)
# pL2 = numpy.linalg.inv(numpy.dot(numpy.dot(pA2.T,numpy.linalg.inv(pC2)),pA2))
# pR2 = numpy.dot(numpy.dot(pA2.T,numpy.linalg.inv(pC2)),pY2)
#
# psig_m1 = numpy.sqrt(pL1[1, 1])
# psig_c1 = numpy.sqrt(pL1[0, 0])
# psig_m2 = numpy.sqrt(pL2[1, 1])
# psig_c2 = numpy.sqrt(pL2[0, 0])
# pe1mc = numpy.dot(pL1, pR1)
# pe2mc = numpy.dot(pL2, pR2)
#
# x = numpy.linspace(-0.02,0.02,100)
# print(me1mc,me2mc)
# print(pe1mc,pe2mc)
#
#
# plt.subplot(121)
# #mean
# labelm1 = 'y=%.6f*x+%.7f'%(me1mc[1],me1mc[0])
# labelp1 = 'y=%.6f*x+%.7f'%(pe1mc[1],pe1mc[0])
# plt.scatter(og1,mean_method[:,0],c='b')
# plt.errorbar(og1,mean_method[:,0],mean_method[:,1],fmt='none',ecolor='b')
# plt.plot(x,me1mc[1]*x+me1mc[0],color='b',label=labelm1)
# #sym-PDF
# plt.scatter(og1,smpdf[:,0],c='r')
# plt.errorbar(og1,smpdf[:,0],smpdf[:,1],fmt='none',ecolor='r')
# plt.plot(x,pe1mc[1]*x+pe1mc[0],color='r',label=labelp1)
# #y=x
# plt.plot([-0.02,0.02],[-0.02,0.02],color='k')
# plt.xlim(-0.015,0.015)
# plt.ylim(-0.015,0.015)
# plt.legend()
# plt.title('g1:200K gals each point')
#
# plt.subplot(122)
# #mean
# labelm2 = 'y=%.6f*x+%.7f'%(me2mc[1],me2mc[0])
# labelp2 = 'y=%.6f*x+%.7f'%(pe2mc[1],pe2mc[0])
# plt.scatter(smpdf[:,-1],mean_method[:,2],c='b')
# plt.errorbar(smpdf[:,-1],mean_method[:,2],mean_method[:,3],fmt='none',ecolor='b')
# plt.plot(x,me2mc[1]*x+me2mc[0],color='b',label = labelm2)
# #sym-PDF
# plt.scatter(smpdf[:,-1],smpdf[:,2],c='r')
# plt.errorbar(smpdf[:,-1],smpdf[:,2],smpdf[:,3],fmt='none',ecolor='r')
# plt.plot(x,pe2mc[1]*x+pe2mc[0],color='r',label= labelp2)
# #y=x
# plt.plot([-0.02,0.02],[-0.02,0.02],color='k')
# #plt.axes().set_aspect('equal', 'datalim')
# plt.xlim(-0.015,0.015)
# plt.ylim(-0.015,0.015)
# plt.legend()
# plt.title('g2:200K gals each point')
# plt.show()




