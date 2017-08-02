import matplotlib.pyplot as plt
import numpy
import time
import pandas
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
mag = Fourier_Quad().divide_stamps(fits.open('E:/gal_chip_12.fits')[0].data,80)

x = []


for i in range(10000):

    pow = Fourier_Quad().pow_spec(mag[i])
    # hlr,pos = Fourier_Quad().get_radius_new(pow,2)
    # signal = 0
    # for a in pos:
    #     signal+= pow[a[0],a[1]]
    #
    noise = numpy.sum(pow*Fourier_Quad().border(1,80))/numpy.sum(Fourier_Quad().border(1,80))
    snr = numpy.sqrt(numpy.mean(pow[39:42,39:42])/noise)
    if snr >500:
        print(i)
    x.append(snr)
    plt.imshow(mag[i])
    title = '%.2f'%snr
    plt.title(title)
    plt.colorbar()
    plt.show()

# plt.scatter(numpy.linspace(0,9999,10000),x)
# plt.show()
# plt.scatter(x,y)
# plt.show()