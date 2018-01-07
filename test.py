from Fourier_Quad import Fourier_Quad
import lsstetc
import numpy
import matplotlib.pyplot as plt
import tool_box
import h5py

data = numpy.zeros((10,10))
f = h5py.File("E:/h5.hdf5","r+")
print(f['/a/data'].value)
# f1 = f.create_group('a')
# print(f1.name)
# f1.create_dataset('data', data=data)
f.close()

# b = [2,5,6,1,5,8]
# a = numpy.array([1,2,6,9,40, 1000])
# plt.plot(a,b)
# plt.xscale("log")
# x1, x2 = plt.xlim()
# plt.plot([x1, x2], [1,1])
# plt.xlim(x1,x2)
#
# plt.show()
# size = 52
# num = 20000
# prop = lsstetc.ETC(band='r', pixel_scale=0.2, stamp_size=size, nvisits=180)
# mags = tool_box.mags_mock(num, 21, 26.5)
#
# noise_sig = prop.sigma_sky
# print(noise_sig)
# fq = Fourier_Quad(size, numpy.random.randint(0, 10000000, 1))
# psf =fq.cre_psf(6, "Moffat")
# plt.imshow(psf)
# plt.show()
# print("PSF:", numpy.sum(psf))
# for k in range(30):
#     p = fq.ran_pos(45, 9)
#     mg1, mg2, mn = 0, 0, 0
#     for i in range(4):
#         p_r = fq.rotate(p, i*numpy.pi/4.)
#         p_s = fq.shear(p_r, -0.01, 0.03)
#         gal_final = fq.convolve_psf(p_s, 6, 100, 'Moffat')
#         a, b, c = fq.shear_est(gal_final, psf, F=False)[0:3]
#         mg1 +=a
#         mg2 +=b
#         mn += c
#     print(mg1/mn, mg2/mn)
#
