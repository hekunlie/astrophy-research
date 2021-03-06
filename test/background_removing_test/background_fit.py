from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
from astropy.io import fits


order = 5
terms = int((order + 1) * (order + 2) / 2)
data_path = "D:/noname/HDtest/"

tag = "bkg1537s"
img_name = "source1%s.fits"%tag

background_og_name = data_path + "%s.fits"%tag
src_og_name = data_path + "source800#1.fits"

background_name = data_path + "source1%s_bg.fits"%tag
background_name_diff = data_path + "source1%s_bg_diff.fits"%tag

src_fit_name = data_path + "source1%s_remove_bg.fits"%tag
src_fit_name_diff = data_path + "source1%s_src_diff.fits"%tag

val_scale = 100000
src_img = fits.open(data_path + img_name)[0].data*val_scale

# nx,ny = 1000,1000
ny, nx = src_img.shape
pix_num = nx*ny
my, mx = numpy.mgrid[0:ny,0:nx]/1000
my_flat = my.flatten()
mx_flat = mx.flatten()
# print(my)
# print(mx)
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

src_img = src_img/val_scale
fit_bg_img = fxy_fit/val_scale
img = Image_Plot()
img.subplots(1,2)
# img.axs[0][0].hist(src_img_flat, 1000)
img.axs[0][1].imshow(fit_bg_img)
img.axs[0][0].imshow(src_img)
img.show_img()

hdu = fits.PrimaryHDU(fit_bg_img)
hdu.writeto(background_name,overwrite=True)

background_og = fits.open(background_og_name)[0].data
diff = background_og - fit_bg_img
hdu = fits.PrimaryHDU(diff)
hdu.writeto(background_name_diff,overwrite=True)

src_fit = src_img - fit_bg_img
hdu = fits.PrimaryHDU(src_fit)
hdu.writeto(src_fit_name,overwrite=True)

src_og = fits.open(src_og_name)[0].data
diff = src_og - src_fit
hdu = fits.PrimaryHDU(diff)
hdu.writeto(src_fit_name_diff,overwrite=True)