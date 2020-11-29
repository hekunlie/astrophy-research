from astropy.io import fits
import numpy
import matplotlib.pyplot as plt


img_ori = fits.open("CSST-PSF-ccd7-w1-1.fits")[0].data
print(img_ori.sum())
size = img_ori.shape[0]
print(img_ori.shape)
print(img_ori.max())
npw = numpy.where(img_ori == img_ori.max())
print(npw)
idx = img_ori >= img_ori.max()/100000
print(img_ori[idx].min())
mask = numpy.zeros_like(img_ori)
mask[idx] = 1
hdu = fits.PrimaryHDU(mask)
hdu.writeto("mask.fits",overwrite=True)
d = 64+16
mask = mask[d:size-d,d:size-d]
# for i in range()
print(mask.shape)
plt.imshow(mask)
plt.show()

hdu = fits.PrimaryHDU(img_ori[d:size-d,d:size-d])
hdu.writeto("CSST_PSF_cut.fits",overwrite=True)