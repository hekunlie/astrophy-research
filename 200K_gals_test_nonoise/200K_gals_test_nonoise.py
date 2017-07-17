from sys import path
#path.append('E:/Github/astrophy-research')
from Fourier_Quad import Fourier_Quad
import numpy
import pandas
import matplotlib.pyplot as plt
<<<<<<< HEAD:point_sources.py


def g_fit(G,N,U,mode,num,left=-0.06,right=0.06):
    g_range = numpy.linspace(left,right,num)
    xisq = [Fourier_Quad().G_bin(G,N,U,g_h,mode,8) for g_h in g_range]
    gg4 = numpy.sum(g_range ** 4)
    gg3 = numpy.sum(g_range ** 3)
    gg2 = numpy.sum(g_range ** 2)
    gg1 = numpy.sum(g_range)
    xigg2 = numpy.sum(xisq * g_range ** 2)
    xigg1 = numpy.sum(xisq * g_range)
    xigg0 = numpy.sum(xisq)
    cov = numpy.linalg.inv(numpy.array([[gg4, gg3, gg2], [gg3, gg2, gg1], [gg2, gg1, num]]))
    paras = numpy.dot(cov, numpy.array([xigg2, xigg1, xigg0]))
    g_sig = numpy.sqrt(1 / 2. / paras[0])
    g_est = -paras[1]/2/paras[0]
    return g_est,g_sig

og1 = numpy.linspace(-0.01,0.01,11)
og2 = numpy.linspace(-0.01,0.01,11)
numpy.random.shuffle(og2)

size = 60
psf_ori = Fourier_Quad().cre_psf(4,size)

mean_method = numpy.zeros((11,6))
smpdf = numpy.zeros((11,6))

for i in range(len(og1)):
    g1 = og1[i]
    g2 = og2[i]
    data_matrix = 0
    p=0
    for m in range(500):
        points= Fourier_Quad().ran_pos(50,size)
        for k in range(4):
            angle = numpy.pi*k/4
            points_rotate = Fourier_Quad().rotate(points,angle)
            points_shear = Fourier_Quad().shear(points_rotate,g1,g2)
            final = Fourier_Quad().convolve_psf(points_shear,4,size)
            G1,G2,N,U = Fourier_Quad().shear_est(final,psf_ori,size,F=False)[0:4]
            ith_row = numpy.array([ G1, G2, N,U,g1,g2])
            if p == 0:
                data_matrix = ith_row
            else:
                data_matrix = numpy.row_stack((data_matrix, ith_row))
            p = 1
=======
from scipy import optimize
import time
import os


mean_method = numpy.zeros((11,6))
smpdf = numpy.zeros((11,6))
# pdf_g = numpy.zeros((11,4))
files = os.listdir('E:/Github/astrophy-research/200K_gals_test_nonoise/')
path = []
for file in files:
    if 'gal_test_' in file:
        path.append(file)

for i in range(len(path)):
    data_matrix = numpy.load(path[i])['arr_0']
>>>>>>> 16cb26ce1feeb597502df492eb0f8d6d6aacc6f2:200K_gals_test_nonoise/200K_gals_test_nonoise.py
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
    g1_est,sig1  = Fourier_Quad().fmin_g(G1,N,U,mode=1,bin_num=8,sample=20000)[0:2]
    g2_est,sig2 = Fourier_Quad().fmin_g(G2, N, U, mode=2, bin_num=8,sample=20000)[0:2]

<<<<<<< HEAD:point_sources.py
    g1_est,sig1  = Fourier_Quad().fmin_g(G1,N,U,mode=1,bin_num=6,left=-0.02,right=0.02,iters=12,sample=None)[0:2]
    g2_est,sig2 = Fourier_Quad().fmin_g(G2, N, U, mode=2, bin_num=6, left=-0.02, right=0.02,iters=12, sample=None)[0:2]
    smpdf[i] = g1_est,sig1,g2_est,sig2,g1,g2
    print(g1,g1_est_mean,g1_est)
    print(g2,g2_est_mean,g2_est)
=======
    #print(g1,g1_est_mean,g2,g2_est_mean)
    print(g1,g1_est,sig1,g2,g2_est,sig2)
    smpdf[i] = g1_est,sig1,g2_est,sig2,g1,g2
    mean_method[i] = g1_est_mean,g1_est_mean_sig,g2_est_mean,g2_est_mean_sig,g1,g2




>>>>>>> 16cb26ce1feeb597502df492eb0f8d6d6aacc6f2:200K_gals_test_nonoise/200K_gals_test_nonoise.py
numpy.savez('cache.npz',mean_method,smpdf)

result = numpy.load('cache.npz')
mean_method = result['arr_0']
smpdf = result['arr_1']
og1= mean_method[:,-2]
og2 = mean_method[:,-1]
mA1 = numpy.column_stack((numpy.ones_like(og1),og1))
mY1  = mean_method[:,0]
mC1  = numpy.diag(mean_method[:,1]**2)
mA2 = numpy.column_stack((numpy.ones_like(smpdf[:,-1]),smpdf[:,-1]))
mY2  = mean_method[:,2]
mC2  = numpy.diag(mean_method[:,3]**2)

mL1 = numpy.linalg.inv(numpy.dot(numpy.dot(mA1.T,numpy.linalg.inv(mC1)),mA1))
mR1 = numpy.dot(numpy.dot(mA1.T,numpy.linalg.inv(mC1)),mY1)
mL2 = numpy.linalg.inv(numpy.dot(numpy.dot(mA2.T,numpy.linalg.inv(mC2)),mA2))
mR2 = numpy.dot(numpy.dot(mA2.T,numpy.linalg.inv(mC2)),mY2)

msig_m1 = numpy.sqrt(mL1[1, 1])
msig_c1 = numpy.sqrt(mL1[0, 0])
msig_m2 = numpy.sqrt(mL2[1, 1])
msig_c2 = numpy.sqrt(mL2[0, 0])
me1mc = numpy.dot(mL1, mR1)
me2mc = numpy.dot(mL2, mR2)

pA1 = numpy.column_stack((numpy.ones_like(og1),og1))
pY1  = smpdf[:,0]
pC1  = numpy.diag(smpdf[:,1]**2)
pA2 = numpy.column_stack((numpy.ones_like(smpdf[:,-1]),smpdf[:,-1]))
pY2  = smpdf[:,2]
pC2  = numpy.diag(smpdf[:,3]**2)

pL1 = numpy.linalg.inv(numpy.dot(numpy.dot(pA1.T,numpy.linalg.inv(pC1)),pA1))
pR1 = numpy.dot(numpy.dot(pA1.T,numpy.linalg.inv(pC1)),pY1)
pL2 = numpy.linalg.inv(numpy.dot(numpy.dot(pA2.T,numpy.linalg.inv(pC2)),pA2))
pR2 = numpy.dot(numpy.dot(pA2.T,numpy.linalg.inv(pC2)),pY2)

psig_m1 = numpy.sqrt(pL1[1, 1])
psig_c1 = numpy.sqrt(pL1[0, 0])
psig_m2 = numpy.sqrt(pL2[1, 1])
psig_c2 = numpy.sqrt(pL2[0, 0])
pe1mc = numpy.dot(pL1, pR1)
pe2mc = numpy.dot(pL2, pR2)

x = numpy.linspace(-0.02,0.02,100)
print(me1mc,me2mc)
print(pe1mc,pe2mc)
# plt.figure(figsize=(20,10))
plt.subplot(121)
#mean
<<<<<<< HEAD:point_sources.py
labelm1 = 'y = %.6f*x+%.9f'%(me1mc[1],me1mc[0])
plt.scatter(og1,mean_method[:,0],c='g')
plt.errorbar(og1,mean_method[:,0],mean_method[:,1],fmt='none',ecolor='g')
plt.plot(x,me1mc[1]*x+me1mc[0],color='g',label=labelm1)
=======
labelm1 = 'y=%.6f*x+%.7f'%(me1mc[1],me1mc[0])
labelp1 = 'y=%.6f*x+%.7f'%(pe1mc[1],pe1mc[0])
plt.scatter(og1,mean_method[:,0],c='b')
plt.errorbar(og1,mean_method[:,0],mean_method[:,1],fmt='none',ecolor='b')
plt.plot(x,me1mc[1]*x+me1mc[0],color='b',label=labelm1)
>>>>>>> 16cb26ce1feeb597502df492eb0f8d6d6aacc6f2:200K_gals_test_nonoise/200K_gals_test_nonoise.py
#sym-PDF
labelp1 = 'y = %.6f*x+%.9f'%(pe1mc[1],pe1mc[0])
plt.scatter(og1,smpdf[:,0],c='r')
plt.errorbar(og1,smpdf[:,0],smpdf[:,1],fmt='none',ecolor='r')
plt.plot(x,pe1mc[1]*x+pe1mc[0],color='r',label=labelp1)
#y=x
plt.plot([-0.02,0.02],[-0.02,0.02],color='k')
#plt.axes().set_aspect('equal', 'datalim')
plt.xlim(-0.015,0.015)
plt.ylim(-0.015,0.015)
plt.legend()
<<<<<<< HEAD:point_sources.py
=======
plt.title('g1:200K gals each point')
plt.subplot(122)
>>>>>>> 16cb26ce1feeb597502df492eb0f8d6d6aacc6f2:200K_gals_test_nonoise/200K_gals_test_nonoise.py

plt.subplot(122)
#mean
<<<<<<< HEAD:point_sources.py
labelm2 = 'y = %.6f*x+%.9f'%(me2mc[1],me2mc[0])
labelp2 = 'y = %.6f*x+%.9f'%(pe2mc[1],pe2mc[0])
plt.scatter(smpdf[:,-1],mean_method[:,2],c='g')
plt.errorbar(smpdf[:,-1],mean_method[:,2],mean_method[:,3],fmt='none',ecolor='g')
plt.plot(x,me2mc[1]*x+me2mc[0],color='g',label = labelm2)
#sym-PDF
plt.scatter(smpdf[:,-1],smpdf[:,2],c='r')
plt.errorbar(smpdf[:,-1],smpdf[:,2],smpdf[:,3],fmt='none',ecolor='r')
plt.plot(x,pe2mc[1]*x+pe2mc[0],color='r',label=labelp2 )
#y=x
plt.plot(x,x,color='k',linestyle='-')
plt.xlim(-0.02,0.02)
plt.ylim(-0.02,0.02)
=======
labelm2 = 'y=%.6f*x+%.7f'%(me2mc[1],me2mc[0])
labelp2 = 'y=%.6f*x+%.7f'%(pe2mc[1],pe2mc[0])
plt.scatter(smpdf[:,-1],mean_method[:,2],c='b')
plt.errorbar(smpdf[:,-1],mean_method[:,2],mean_method[:,3],fmt='none',ecolor='b')
plt.plot(x,me2mc[1]*x+me2mc[0],color='b',label = labelm2)
#sym-PDF
plt.scatter(smpdf[:,-1],smpdf[:,2],c='r')
plt.errorbar(smpdf[:,-1],smpdf[:,2],smpdf[:,3],fmt='none',ecolor='r')
plt.plot(x,pe2mc[1]*x+pe2mc[0],color='r',label= labelp2)
#y=x
plt.plot([-0.02,0.02],[-0.02,0.02],color='k')
#plt.axes().set_aspect('equal', 'datalim')
plt.xlim(-0.015,0.015)
plt.ylim(-0.015,0.015)
>>>>>>> 16cb26ce1feeb597502df492eb0f8d6d6aacc6f2:200K_gals_test_nonoise/200K_gals_test_nonoise.py

plt.legend()
plt.title('g2:200K gals each point')
plt.show()




