from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
from astropy.io import fits


order = 5
terms = int((order + 1) * (order + 2) / 2)
val_scale = 100000
src_img = fits.open("D:/ll800+fg.fits")[0].data*val_scale
# nx,ny = 1000,1000
ny, nx = src_img.shape
pix_num = nx*ny
my, mx = numpy.mgrid[0:ny,0:nx]/1000
my_flat = my.flatten()
mx_flat = mx.flatten()
print(my)
print(mx)
print(src_img.max(), src_img.min())
# noise = numpy.random.normal(0,10,(ny,nx))
# src_img = 10 + mx*2 + my*0.001 + mx*my*0.001 + mx**2*0.00005 + noise

src_img_flat = src_img.flatten()
seq = numpy.arange(0,pix_num)
pix_lw = src_img_flat[int(pix_num/50)]
pix_up = src_img_flat[int(pix_num*49/50)]
img = Image_Plot()
img.subplots(1,1)
img.axs[0][0].hist(src_img_flat, 1000)
# img.axs[0][0].imshow(src_img)
img.show_img()
fit_para = tool_box.fit_background(src_img, 300000, "flat", pix_lw, pix_up, my_flat, mx_flat, seq, order=order)[0]
print(fit_para)
fxy_fit = 0
pows = [(int(i-j), int(j)) for i in range(order+1) for j in range(i+1)]
for i in range(terms):
    fxy_fit += fit_para[i,0]*mx**(pows[i][0])*my**(pows[i][1])

img = Image_Plot()
img.subplots(1,2)
# img.axs[0][0].hist(src_img_flat, 1000)
img.axs[0][1].imshow(fxy_fit/val_scale)
img.axs[0][0].imshow(src_img/val_scale)
img.show_img()
hdu = fits.PrimaryHDU(fxy_fit/val_scale)
hdu.writeto("D:/background.fits",overwrite=True)

hdu = fits.PrimaryHDU(src_img/val_scale - fxy_fit/val_scale)
hdu.writeto("D:/source.fits",overwrite=True)