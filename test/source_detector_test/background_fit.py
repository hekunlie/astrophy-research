from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
from astropy.io import fits

expo = fits.open("D:/frame-1.fits")[0].data
print(expo.shape)
expo_pixel = numpy.sort(expo.flatten())
sigma = expo_pixel.std()
mean = expo_pixel.mean()
pixel_num = len(expo_pixel)
idx1 = expo_pixel >= -0.5#expo_pixel[int(pixel_num/20)]#mean - sigma*2
idx2 = expo_pixel <= 0.1#expo_pixel[int(pixel_num/20*19)]#mean + sigma*10
idx = idx1 & idx2
print(expo_pixel[int(pixel_num/20)], expo_pixel[int(pixel_num*19/20)])
print(numpy.std(expo_pixel[idx]))
# idx = expo_pixel < 0.15
bin_num = 50
num, bins = numpy.histogram(expo_pixel[idx], bin_num)[:2]
x = (bins[:bin_num] + bins[1:])/2


print(num)
img = Image_Plot()
img.subplots(1,2)
img.axs[0][0].imshow(expo)
img.axs[0][1].hist(expo_pixel[idx],bin_num)
img.show_img()
