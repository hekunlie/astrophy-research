from Fourier_Quad import Fourier_Quad
import lsstetc
import numpy
from astropy.io import fits
import matplotlib.pyplot as plt
import tool_box
import h5py
import copy
import time
from scipy import optimize

size = 52
fq = Fourier_Quad(size, numpy.random.randint(0, 10000000, 1))
# img = fits.open("E:/gal_chip_0039.fits")[0].data
# c_img = copy.copy(img)
#
# pixels = numpy.sort(numpy.ndarray.flatten(img))
# #pixels = pixels[pixels>-2000]
# bins, num = plt.hist(pixels[0:int(0.99*len(pixels))], 100)[0:2]
# #plt.show()
# x_num = len(num)
# print(len(bins))
# print(len(num))
# xx = numpy.linspace(-1000, 1000, 200)
#
# def f(x, a, b, c):
#     return c*numpy.exp(-(x-a)**2/2/b**2)
#
#
# p, pc = optimize.curve_fit(f, num[:x_num-1], bins, method='lm')
# plt.close()
# plt.plot(xx, f(xx, p[0], p[1], p[2]))
# plt.show()
# print(p)
# idx = c_img < 2*p[1]
# c_img[idx] = 0
# idx = c_img > 0
# c_img[idx] = 1
#
# img_2 = copy.copy(c_img)
# y, x = c_img.shape
# print(y,x)
# t1 = time.clock()
# gals = tool_box.source_detector(c_img, y, x)
# t2 = time.clock()
#
# print(len(gals),t2-t1)
#
# fig = numpy.zeros_like(c_img)
# for gal in gals:
#     for cor in gal:
#         fig[cor[0], cor[1]] = img_2[cor[0], cor[1]]
#
# plt.subplot(131)
# plt.imshow(img)
# plt.subplot(132)
# plt.imshow(img_2)
# plt.subplot(133)
# plt.imshow(fig)
# plt.show()
#
# flux = []
# m = 0
# t = -1
# gal = gals[100]
# print(gal)
# for i, cor in enumerate(gal):
#     if img[cor[0], cor[1]] > m:
#         m = img[cor[0], cor[1]]
#         t = i
#         print(m,t)
# print(t,m)
# y,x = gal[t]
# plt.imshow(img[y-25:y+25, x-25:x+25])
# plt.show()



# psf =fq.cre_psf(4, "Moffat")
# plt.imshow(psf)
# plt.show()
# print("PSF:", numpy.sum(psf))
# pp = fq.pow_spec(psf)
# plt.imshow(pp)
# plt.show()
# r1 = fq.get_radius_new(pp, 2)[0]
# my, mx = numpy.mgrid[0:size, 0:size] - size/2
# wb1 = numpy.exp(-(mx**2+my**2)/r1**2)
#
# r2 = fq.get_radius_new(pp, 2.7)[0]/2
# wb2 = numpy.exp(-2*(mx**2+my**2)/r2**2)
# plt.subplot(221)
# plt.imshow(pp)
# plt.subplot(222)
# plt.imshow(wb1)
# plt.subplot(223)
# plt.imshow(wb2)
# plt.show()

# res = numpy.zeros((num,4))
for k in range(10):
    p = fq.ran_pos(45, 9)
    p_s = fq.shear(p, -0.01, 0.03)
    noise = fq.draw_noise(0,20)
    gal_final = fq.convolve_psf(p_s, 4, 1000000, 'Moffat')# + fq.draw_noise(0, 40)
    print(numpy.sum(gal_final))
#    # plt.imshow(gal_final)
#    # plt.show()
#     a, b, c, d= fq.shear_est(gal_final, psf, noise, F=False)[0:4]
#     res[k,0] = a
#     res[k,1] = b
#     res[k,2] = c
#     res[k,3] = d
# g, gh = fq.fmin_g(g=res[:,0], n=res[:,2], u=res[:,3], mode=1, bin_num=8)
# print(g, gh)
# print(numpy.sqrt(num)*gh)

