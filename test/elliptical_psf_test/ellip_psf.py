import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot


h5f = h5py.File("./img.hdf5","r")
epsf_cpp = h5f["/epsf"][()]
gal_epsf_cpp = h5f["/gal_epsf"][()]
cpsf_cpp = h5f["/cpsf"][()]
gal_cpsf_cpp = h5f["/gal_cpsf"][()]
pts = h5f["/pts"][()]
# h5f.close()

psf_scale = 4
size = 48
psf_e = 0.6
theta = 0.5*numpy.pi
fq = Fourier_Quad(size,112)
pts_num = pts.shape[1]
# print(pts_num)
# print(pts)
# print("C-PSF")
# cpsf = fq.cre_psf(psf_scale,ellip_theta=(0, 0))
# gal = fq.convolve_psf(pts, psf_scale, 1*pts_num, ellip_theta=(0, 0))
# x,y = fq.mx.flatten(), fq.my.flatten()
# z1,z2 = cpsf.flatten(),cpsf_cpp.flatten()
# print(z1.sum(), z2.sum())
# gal_z1, gal_z2 = gal.flatten(), gal_cpsf_cpp.flatten()
# diff_psf = z1 - z2
# diff_gal = gal_z1 - gal_z2
#
# print("PSF diff: ", diff_psf.max(), diff_psf.min())
# print("Gal diff: ", diff_gal.max(), diff_gal.min())
#
# img = Image_Plot()
# img.subplots(2, 3)
# img.scatter_pts(0,0,x,y,z1,pts_size=15,color_map="jet",marker="o")
# img.scatter_pts(0,1,x,y,z2,pts_size=15,color_map="jet",marker="o")
# img.scatter_pts(0,2,x,y,diff_psf,pts_size=15,color_map="jet",marker="o")
# img.scatter_pts(1,0,x,y,gal_z1,pts_size=15,color_map="jet",marker="o")
# img.scatter_pts(1,1,x,y,gal_z2,pts_size=15,color_map="jet",marker="o")
# img.scatter_pts(1,2,x,y,diff_gal,pts_size=15,color_map="jet",marker="o")
# # img.scatter_pts(2,i,x,y,diff,pts_size=25,color_map="bwr")
# img.show_img()



print("E-PSF")
gal_imgs = numpy.zeros((size, 4*size))
# gal_imgs_cpp = numpy.zeros((size, 4*size))
psf_imgs = numpy.zeros((size, 4*size))
# psf_imgs_cpp = numpy.zeros((size, 4*size))

for i in range(4):
    psf = fq.cre_psf(psf_scale,ellip_theta=(psf_e, numpy.pi/4*i + theta))
    gal = fq.convolve_psf(pts, psf_scale, 1, ellip_theta=(psf_e, numpy.pi/4*i + theta))
    gal_imgs[:,int(i*size):int((i+1)*size)] = gal
    psf_imgs[:,int(i*size):int((i+1)*size)] = psf

diff_psf = psf_imgs - epsf_cpp
print(diff_psf.min(), diff_psf.max())
img = Image_Plot(fig_x=4,fig_y=3)
img.subplots(3, 1)
fig = img.axs[0][0].imshow(psf_imgs)
img.figure.colorbar(fig,ax=img.axs[0][0])
fig = img.axs[1][0].imshow(epsf_cpp)
img.figure.colorbar(fig,ax=img.axs[1][0])
fig = img.axs[2][0].imshow(diff_psf)
img.figure.colorbar(fig,ax=img.axs[2][0])
img.show_img()

diff_gal = gal_imgs - gal_epsf_cpp
print(diff_gal.min(), diff_gal.max())

img = Image_Plot(fig_x=4,fig_y=3)
img.subplots(3, 1)
fig = img.axs[0][0].imshow(gal_imgs)
img.figure.colorbar(fig,ax=img.axs[0][0])
fig = img.axs[1][0].imshow(gal_epsf_cpp)
img.figure.colorbar(fig,ax=img.axs[1][0])
fig = img.axs[2][0].imshow(diff_gal)
img.figure.colorbar(fig,ax=img.axs[2][0])
img.show_img()
