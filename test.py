from mpi4py import MPI
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import numpy
from Fourier_Quad import Fourier_Quad
import time
from astropy.io import fits
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


fg1 = numpy.linspace(-0.02, 0.02, 8)
fg2 = -numpy.linspace(-0.02, 0.02, 8)
scale = 4
size = 52
num = 10000
chip_num = 10
fq = Fourier_Quad(size, (rank+1)*2441452)

psf = fq.cre_psf(scale, "Moffat")
psfp = fq.pow_spec(psf)
g1, g2 = fg1[rank], fg2[rank]
estor = numpy.zeros((num*chip_num, 5))
data_path = "/lmc/test/data/data_%d.txt" %rank

for i in range(chip_num):
    t1 = time.clock()
    chip_path = "/lmc/test/%d_gal_chip_%d.fits" %(rank, i)
    img = []
    for j in range(num):
        noise1 = fq.noise(0, 5)
        p = fq.ran_pos(45, 10, (g1, g2))[1]
        gal = fq.convolve_psf(p, scale, flux=10, psf="Moffat") + noise1
        noise2 = fq.noise(0, 5)
        mg1, mg2, mn = fq.shear_est(gal, psfp, noise2, F=True)[0:3]
        estor[j+i*num, 0] = mg1
        estor[j+i*num, 1] = mg2
        estor[j+i*num, 2] = mn
        estor[j+i*num, 3] = g1
        estor[j+i*num, 4] = g2
        if rank == 0 and i == 0:
            img.append(gal)

    if rank == 0 and i == 0:
        chip = fq.image_stack(img, 100)
        hdu = fits.PrimaryHDU(chip)
        hdu.writeto(chip_path, overwrite=True)
    t2 = time.clock()

numpy.savetxt(data_path, estor)
estg1 = numpy.mean(estor[:, 0]) / numpy.mean(estor[:, 2])
estdg1 = numpy.sqrt(numpy.mean(estor[:, 0] ** 2) / numpy.mean(estor[:, 2]) ** 2 / len(estor[:, 2]))
estg2 = numpy.mean(estor[:, 1]) / numpy.mean(estor[:, 2])
estdg2 = numpy.sqrt(numpy.mean(estor[:, 1] ** 2) / numpy.mean(estor[:, 2]) ** 2 / len(estor[:, 2]))
print(rank, estg1, estdg1, estg2, estdg2)

if rank == 0:
    est_g = numpy.zeros((4, 8))
    for i in range(8):
        data_path = "/lmc/test/data/data_%d.txt" % i
        data = numpy.loadtxt(data_path)
        mg1 = data[:, 0]
        mg2 = data[:, 1]
        mn = data[:, 2]

        est_g[0,i] = numpy.mean(mg1)/numpy.mean(mn)
        est_g[1,i] = numpy.sqrt(numpy.mean(mg1**2) / numpy.mean(mn)**2/len(mg1))

        est_g[2,i] = numpy.mean(mg2)/numpy.mean(mn)
        est_g[3,i] = numpy.sqrt(numpy.mean(mg2**2) / numpy.mean(mn)**2/len(mg2))
        print(i, len(mg1), fg1[i], est_g[0,i], est_g[1,i])
        print(i, len(mg2), fg2[i], est_g[2,i], est_g[3,i])

    A1 = numpy.column_stack((numpy.ones_like(fg1.T), fg1.T))
    Y1 = est_g[0].T
    C1 = numpy.diag(est_g[1]**2)

    A2 = numpy.column_stack((numpy.ones_like(fg2.T), fg2.T))
    Y2 = est_g[2].T
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
    plt.errorbar(fg1, est_g[0], est_g[1])
    plt.plot(fg1, fg1)
    plt.subplot(122)
    plt.errorbar(fg2, est_g[2], est_g[3])
    plt.plot(fg2, fg2)
    plt.savefig("/home/hklee/test.png")