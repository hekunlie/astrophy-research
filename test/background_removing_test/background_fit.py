from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
from astropy.io import fits


order = 3
terms = int((order + 1) * (order + 2) / 2)

src_img = fits.open("D:/ll800+fg.fits")[0].data*100000

ny, nx = src_img.shape
pix_num = nx*ny
my, mx = numpy.mgrid[0:ny,0:nx]
my_flat = my.flatten()
mx_flat = mx.flatten()
print(my)
print(mx)
noise = numpy.random.normal(0,10,(ny,nx))
src_img = 10 + mx*2 + my*0.001 + mx*my*0.001 + mx**2*0.00005 + noise

src_img_flat = src_img.flatten()
seq = numpy.arange(0,pix_num)
pix_lw = src_img_flat[int(pix_num/5)]
pix_up = src_img_flat[int(pix_num*4/5)]
img = Image_Plot()
img.subplots(1,1)
img.axs[0][0].hist(src_img_flat, 1000)
img.show_img()
fit_para = tool_box.fit_background(src_img, 100000, "flat", pix_lw, pix_up, my_flat, mx_flat, seq, order=order)[0]
print(fit_para)
fxy_fit = 0
pows = [(int(i-j), int(j)) for i in range(order+1) for j in range(i+1)]
for i in range(terms):
    fxy_fit += fit_para[i,0]*mx**(pows[i][0])*my**(pows[i][1])

hdu = fits.PrimaryHDU(fxy_fit)
hdu.writeto("D:/background.fits",overwrite=True)

