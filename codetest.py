from sys import path
path.append('/home/hklee/work/fourier_quad/')
import numpy
from Fourier_Quad import Fourier_Quad
import time
import matplotlib.pyplot as plt


fg1 = numpy.linspace(-0.02, 0.02, 8)
fg2 = -numpy.linspace(-0.02, 0.02, 8)
scale = 4
size = 56
num = 10000
psf = Fourier_Quad().cre_psf(scale, size, "Moffat")
psfp = Fourier_Quad().pow_spec(psf)
est_g = numpy.zeros((4, 8))
for i in range(8):
    g1, g2 = fg1[i], fg2[i]
    estor = numpy.zeros((3, num))
    t1 = time.clock()
    for j in range(num):
        noise1 = numpy.random.normal(0, 10, size*size).reshape(size, size)
        noise2 = numpy.random.normal(0, 10, size * size).reshape(size, size)
        p = Fourier_Quad().ran_pos(45, 8, (g1, g2))[1]
        gal = Fourier_Quad().convolve_psf(p, scale, size, flux=10, psf="Moffat") + noise1
        mg1, mg2, mn = Fourier_Quad().shear_est(gal, psfp, size, noise2, F=True)[0:3]
        estor[0, j] = mg1
        estor[1, j] = mg2
        estor[2, j] = mn
    est_g[0, i] = numpy.mean(estor[0]) / numpy.mean(estor[2])
    est_g[1, i] = numpy.mean(estor[1]) / numpy.mean(estor[2])
    est_g[2, i] = numpy.sqrt(numpy.mean(estor[0]**2) / numpy.mean(estor[2])**2/num)
    est_g[3, i] = numpy.sqrt(numpy.mean(estor[1] ** 2) / (numpy.mean(estor[2])) ** 2 / num)
    t2 = time.clock()
    print(t2-t1)

A1 = numpy.column_stack((numpy.ones_like(fg1.T), fg1.T))
Y1 = est_g[0].T
C1 = numpy.diag(est_g[2]**2)

A2 = numpy.column_stack((numpy.ones_like(fg2.T), fg2.T))
Y2 = est_g[1].T
C2 = numpy.diag(est_g[3]**2)

L1 = numpy.linalg.inv(numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), A1))
R1 = numpy.dot(numpy.dot(A1.T, numpy.linalg.inv(C1)), Y1)
L2 = numpy.linalg.inv(numpy.dot(numpy.dot(A2.T, numpy.linalg.inv(C2)), A2))
R2 = numpy.dot(numpy.dot(A2.T, numpy.linalg.inv(C2)), Y2)

sig_m1 = numpy.sqrt(L1[1, 1])
sig_c1 = numpy.sqrt(L1[0, 0])
sig_m2 = numpy.sqrt(L2[1, 1])
sig_c2 = numpy.sqrt(L2[0, 0])
e1mc = numpy.dot(L1, R1)
e2mc = numpy.dot(L2, R2)

print(e1mc[1] - 1, sig_m1,  e1mc[0], sig_c1)
print(e2mc[1] - 1, sig_m2,  e2mc[0], sig_c2)
plt.subplot(121)
plt.errorbar(fg1, est_g[0], est_g[2])
plt.plot(fg1, fg1)
plt.subplot(122)
plt.errorbar(fg2, est_g[1], est_g[3])
plt.plot(fg2, fg2)
plt.show()
