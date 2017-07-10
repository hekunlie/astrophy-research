from Fourier_Quad import Fourier_Quad
import numpy
import pandas
import matplotlib.pyplot as plt

og1 = numpy.linspace(-0.01,0.01,21)
og2 = numpy.linspace(-0.01,0.01,21)
print(og1.shape)
numpy.random.shuffle(og2)

size = 60
psf_ori = Fourier_Quad().cre_psf(4,size)

mean_method = numpy.zeros((21,6))
smpdf = numpy.zeros((21,6))
for i in range(len(og1)):
    g1 = og1[i]
    g2 = og2[i]
    data_matrix = 0
    p=0
    for m in range(2000):
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
    G1 = data_matrix[:,0]
    G2 = data_matrix[:,1]
    N = data_matrix[:,2]
    U = data_matrix[:,3]

    g1_est_mean = numpy.sum(G1)/numpy.sum(N)

    g1_est_mean_sig = numpy.std(G1/N)/numpy.sqrt(len(N))
    g2_est_mean = numpy.sum(G2) / numpy.sum(N)
    g2_est_mean_sig = numpy.std(G2 / N) / numpy.sqrt(len(N))
    mean_method[i] = g1_est_mean,g1_est_mean_sig,g2_est_mean,g2_est_mean_sig,g1,g2

    g1_est,sig1  = Fourier_Quad().fmin_g(G1,N,U,mode=1,bin_num=6,left=-0.02,right=0.02,iters=5,sample=None)[0:2]
    g2_est,sig2 = Fourier_Quad().fmin_g(G2, N, U, mode=2, bin_num=6, left=-0.02, right=0.02,iters=5, sample=None)[0:2]
    smpdf[i] = g1_est,sig1,g2_est,sig2,g1,g2
    # print(g1,g1_est_mean,g1_est)
    # print(g2,g2_est_mean,g2_est)
numpy.savez('cache.npz',mean_method,smpdf)

result = numpy.load('cache.npz')
mean_method = result['arr_0']
smpdf = result['arr_1']

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

psig_m1 = numpy.sqrt(mL1[1, 1])
psig_c1 = numpy.sqrt(mL1[0, 0])
psig_m2 = numpy.sqrt(mL2[1, 1])
psig_c2 = numpy.sqrt(mL2[0, 0])
pe1mc = numpy.dot(mL1, mR1)
pe2mc = numpy.dot(mL2, mR2)

x = numpy.linspace(-0.02,0.02,100)
print(me1mc,me2mc)
print(pe1mc,pe2mc)
# plt.figure(figsize=(20,10))
plt.subplot(121)
#mean
plt.scatter(og1,mean_method[:,0],c='g')
plt.errorbar(og1,mean_method[:,0],mean_method[:,1],fmt='none',ecolor='g')
plt.plot(x,me1mc[1]*x+me1mc[0],color='g',label='mean')
#sym-PDF
plt.scatter(og1,smpdf[:,0],c='r')
plt.errorbar(og1,smpdf[:,0],smpdf[:,1],fmt='none',ecolor='r')
plt.plot(x,pe1mc[1]*x+pe1mc[0],color='r',label='sym-PDF')
#y=x
plt.plot(x,x,color='k')
#plt.axes().set_aspect('equal', 'datalim')
plt.xlim(-0.02,0.02)
plt.ylim(-0.02,0.02)
plt.legend()
plt.subplot(122)

#mean
plt.scatter(smpdf[:,-1],mean_method[:,2],c='g')
plt.errorbar(smpdf[:,-1],mean_method[:,2],mean_method[:,3],fmt='none',ecolor='g')
plt.plot(x,me2mc[1]*x+me2mc[0],color='g',label = 'mean')
#sym-PDF
plt.scatter(smpdf[:,-1],smpdf[:,2],c='r')
plt.errorbar(smpdf[:,-1],smpdf[:,2],smpdf[:,3],fmt='none',ecolor='r')
plt.plot(x,pe2mc[1]*x+pe2mc[0],color='r',label= 'sym-PDF',linestyle='-.')
#y=x
plt.plot(x,x,color='k',linestyle='-')
#plt.axes().set_aspect('equal', 'datalim')
plt.xlim(-0.02,0.02)
plt.ylim(-0.02,0.02)

plt.legend()
plt.show()




# df = pandas.DataFrame(data_matrix )
# df.to_excel('E:/data.xlsx')

