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
size = 64
num_p = 45
gal_num = 4
seed = 71#numpy.random.randint(0,100,1)[0]*rank
fq = Fourier_Quad(size,seed)
print("Seed: %d"%seed)
g1 = [0.012, 0.02, -0.01, -0.021, 0]
g2 = [-0.01, -0.015, 0, -0.01, -0.032]

psf_r = numpy.linspace(3, 4.5, 5)

psf = fq.cre_psf(psf_r[2], 600)
# hdu = fits.PrimaryHDU(psf)
# hdu.writeto("E:/psf.fits",overwrite=True)

res = numpy.zeros((4, gal_num))
ppsf0 = fq.pow_spec(psf)
ppsf = tool_box.smooth(ppsf0, size)
plt.subplot(121)
plt.imshow(ppsf)
plt.subplot(122)
plt.imshow(ppsf0)
plt.show()
p = fq.ran_pos(num_p, 8)

sig = 5
for j in range(gal_num):
    t1 = time.clock()
    p_r = fq.rotate(p,j*numpy.pi/4)
    p_s = fq.shear(p_r,g1[1],g2[1])
    gal = fq.convolve_psf(p_s, psf_r[2], 2000)+fq.draw_noise(0, sig)
    noise = fq.draw_noise(0, sig)
    t2 = time.clock()
    plt.imshow(gal)
    plt.show()
    print(j)
    res[0,j], res[1,j], res[2,j], res[3,j] = fq.shear_est(gal, ppsf, noise,F=True)[0:4]
    t3 = time.clock()

print(g1[1], numpy.mean(res[0])/numpy.mean(res[2]))
print(g2[1], numpy.mean(res[1])/numpy.mean(res[2]))
# est_g1, g1_sig = fq.fmin_g(res[0], res[2], res[3], 1, 2)
# est_g2, g2_sig = fq.fmin_g(res[1], res[2], res[3], 2, 2)
#
# if est_g1-2*g1_sig <= g1[rank] <= est_g1+2*g1_sig:
#     t1 = "Consistent with input"
# else:
#     t1 = "out of 2-sigma bound"
# print("rank:%2d, g1: %.5f, est_g1: %.5f (%.5f), %.3f, %s"%(rank,g1[rank], est_g1, g1_sig, numpy.sqrt(gal_num)*g1_sig, t1))
#
# if est_g2-2*g2_sig <= g2[rank] <= est_g2+2*g2_sig:
#     t1 = "Consistent with input"
# else:
#     t1 = "out of 2-sigma bound"
# print("rank:%2d, g2: %.5f, est_g2: %.5f (%.5f), %.3f, %s"%(rank, g2[rank], est_g2, g2_sig, numpy.sqrt(gal_num)*g2_sig, t1))
#




# check the fitting of c++ code


# def f(x, a, b, c, d, e, f):
#     return a * x[0] ** 2 + b * x[0] * x[1] + c * x[1] ** 2 + d * x[0] + e * x[1] + f
#
#
# # gal = fits.open("E:/gal.fits")[0].data
# # gal_pow = numpy.log10(fq.pow_spec(gal))
# # plt.imshow(gal_pow)
# # plt.show()
# # hdu = fits.PrimaryHDU(gal_pow)
# # hdu.writeto("E:/gal_pow.fits",overwrite=True)
# size = 90
# gal_pow = fits.open("E:/gal_pow.fits")[0].data
# my1, mx1 = numpy.mgrid[-2:3,-2:3]
# x, y = mx1.reshape((1, 25)), my1.reshape((1, 25))
# # xx = numpy.delete(x, [0, 4, 20, 24])
# # yy = numpy.delete(y, [0, 4, 20, 24])
#
#
# cen = int((size*size + size)/2)
# fit_img = numpy.zeros_like(gal_pow)
# a = int(size/2)
# for i in range(size):
#     for j in range(size):
#         arr = numpy.zeros((size, size))
#         arr[int(size/2),int(size/2)] = 0.5
#         pos = []
#         tag = 0
#         pk = 0
#         z = []
#         x_d = []
#         for m in range(-2, 3):
#             p = (i + m + size) % size
#             for n in range(-2,3):
#                 q = (j + n + size) % size
#
#                 if tag not in [0,4,20,24]:#abs(m) != 2 or abs(n) != 2:
#                     if p*size+q != cen:
#                         pos.append((p,q))
#                         z.append(gal_pow[p,q])
#                         x_d.append((x[0,tag],y[0,tag]))
#                     else:
#                         pk = tag
#                     #     z.append(0)
#                     #     x_d[0, tag] = 0
#                     #     x_d[1, tag] = 0
#                 # else:
#                 #     z.append(0)
#                 #     x_d[0, tag] = 0
#                 #     x_d[1, tag] = 0
#                 tag+=1
#
#         x_d = numpy.array(x_d).T
#         # print("x:",x_d[0])
#         # print("y:",x_d[1])
#         # s = ""
#
#         # print(pk,len(pos))
#         # for xy in pos:
#         #     arr[xy[0], xy[1]] = 1
#         # plt.imshow(arr[a-10:a+10,a-10:a+10])
#         # plt.grid(which='minor')
#         # plt.show()
#         a1,a2 = curve_fit(f,x_d,z)
#         #print(a1)
#         #print(a2)
#         if i==a and j==a:
#             s = " "
#             for k in range(len(pos)):
#                 s += "%.5f\n " % z[k]
#             print(pk, "z: %s" % s)
#             print(a1[5])
#             print(gal_pow[a, a])
#         fit_img[i,j] = a1[5]
# #fimg = tool_box.smooth(gal_pow, size)
# c_fit = fits.open("E:/galpf.fits")[0].data
#
# plt.subplot(221)
# plt.imshow(gal_pow)
# plt.subplot(222)
# plt.imshow(fit_img)
# plt.subplot(223)
# plt.imshow(c_fit)
# plt.subplot(224)
# plt.imshow(c_fit-fit_img)
# plt.colorbar()
# plt.show()

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
