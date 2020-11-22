from astropy.io import fits
import numpy

img_ori = fits.open("CSST-PSF-ccd7-w1-1.fits")[0].data
print(img_ori.sum())
size = img_ori.shape[0]
print(img_ori.shape)

npw = numpy.where(img_ori == img_ori.max())
print(npw)
idx = img_ori >= img_ori.max()/100000

mask = numpy.zeros_like(img_ori)
mask[idx] = 1
hdu = fits.PrimaryHDU(mask)
hdu.writeto("mask.fits",overwrite=True)
# for i in range()
