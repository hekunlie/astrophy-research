import numpy
from sys import path
path.append('/home/hklee/work/fourier_quad/')
import time
from Fourier_Quad import Fourier_Quad
from sys import argv
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
from scipy.optimize import curve_fit
import tool_box
# from mpi4py import MPI
#
#
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# cpus = comm.Get_size()
rank = 1
size = 90
num_p = 45
gal_num = 1
fq = Fourier_Quad(size,numpy.random.randint(0,100,1)[0]*rank)

g1 = [0.012, 0.02, -0.01, -0.021, 0]
g2 = [-0.01, -0.015, 0, -0.01, -0.032]

psf_r = numpy.linspace(3, 4.5, 5)

psf = fq.cre_psf(psf_r[rank], 600)
hdu = fits.PrimaryHDU(psf)
hdu.writeto("E:/psf.fits",overwrite=True)

res = numpy.zeros((4,gal_num))
ppsf = fq.pow_spec(psf)
for j in range(gal_num):
    p = fq.ran_pos(num_p, 8, (g1[rank], g2[rank]))[1]
    gal = fq.convolve_psf(p, psf_r[rank], 2000)+fq.draw_noise(0,100)
    hdu = fits.PrimaryHDU(gal)
    hdu.writeto("E:/gal.fits", overwrite=True)
    galp= fq.pow_spec(gal)
    plt.subplot(121)
    plt.imshow(gal)
    plt.subplot(122)
    plt.imshow(galp)
    plt.show()
    res[0,j], res[1,j], res[2,j], res[3,j] = fq.shear_est(gal, psf)[0:4]

est_g1, g1_sig = fq.fmin_g(res[0], res[2], res[3], 1, 2)
est_g2, g2_sig = fq.fmin_g(res[1], res[2], res[3], 2, 2)

if est_g1-2*g1_sig <= g1[rank] <= est_g1+2*g1_sig:
    t1 = "Consistent with input"
else:
    t1 = "out of 2-sigma bound"
print("rank:%2d, g1: %.5f, est_g1: %.5f (%.5f), %.3f, %s"%(rank,g1[rank], est_g1, g1_sig, numpy.sqrt(gal_num)*g1_sig, t1))

if est_g2-2*g2_sig <= g2[rank] <= est_g2+2*g2_sig:
    t1 = "Consistent with input"
else:
    t1 = "out of 2-sigma bound"
print("rank:%2d, g2: %.5f, est_g2: %.5f (%.5f), %.3f, %s"%(rank, g2[rank], est_g2, g2_sig, numpy.sqrt(gal_num)*g2_sig, t1))



#
# f = numpy.zeros((1000,1))
# for i in range(1000):
#     n1 = fq.draw_noise(0,380)
#     n = fq.draw_noise(0,380)
#     pos = fq.ran_pos(num_p, 8)
#     gal = fq.convolve_psf(pos, 4, 2000) + n
#     p1 = fq.pow_spec(gal)
#     p2 = fq.pow_spec(n1)
#     p = p1 - p2
#     # plt.subplot(221)
#     # plt.imshow(p1)
#     # plt.subplot(222)
#     # plt.imshow(p2)
#     # plt.subplot(223)
#     plt.imshow(gal)
#     # plt.subplot(224)
#     plt.show()
#     rim = fq.border(1)
#     noise = numpy.sum(rim*(p2-p1))/numpy.sum(rim)
#     if noise<=0 or p[32,32]<=0:
#         print(noise,p[32,32])
#         plt.subplot(131)
#         plt.imshow(gal)
#         plt.subplot(132)
#         plt.imshow(p1)
#         plt.subplot(133)
#         plt.imshow(p)
#         plt.show()
#         plt.close()
#     f[i] = numpy.sqrt(p[32,32]/noise)
# print(f.shape)
# plt.hist(f[:,0],40)
# plt.show()

# check the fitting of c++ code

# path1 = "F:/fit_checking/gal_chip_0000_fit_pow.fits"
# path2 = "F:/fit_checking/gal_chip_0000_pow.fits"
#
# fit_pow = fq.segment(fits.open(path1)[0].data)
# ori_pow = fq.segment(fits.open(path2)[0].data)
#
# for k in range(10000):
#     p10 = ori_pow[k]
#     fit_img = tool_box.ps_fit(p10,size)
#     dif = (fit_pow[k]-fit_img)#*rim
#     plt.subplot(221)
#     plt.imshow(ori_pow[k])
#     plt.subplot(222)
#     plt.imshow(fit_img)
#     plt.subplot(223)
#     plt.imshow(fit_pow[k])
#     plt.subplot(224)
#     plt.imshow(dif)
#     plt.show()
