import numpy
from sys import path
path.append("E:/Github/astrophy-research/my_lib")
from Fourier_Quad import Fourier_Quad
import tool_box
from astropy.io import fits
import matplotlib.pyplot as plt
import h5py
from scipy.optimize import curve_fit

size = 64
f = h5py.File("E:/data.hdf5")
data = f["/data"].value
f.close()
print(data.shape)
flux = data[:,5]
fq = Fourier_Quad(size, 123)
img = fits.open("E:/gal_chip_0010.fits")[0].data
gals = fq.segment(img)
flux_ = []
rim = fq.border(1)
num = rim.sum()

my, mx = numpy.mgrid[-2:3, -2:3]
y, x = numpy.delete(my, [0, 4, 12, 20, 24]).reshape((1,20)),numpy.delete(mx, [0,4, 12, 20, 24]).reshape((1,20))

x_d = numpy.row_stack((y, x))
print(y)
print(x)
for i in range(10000):
    flux_.append(numpy.abs(numpy.sum(gals[i])))
    continue
    gal_p = fq.pow_spec(gals[i])
    sub_gal = numpy.log10(gal_p[int(size/2-2):int(size/2+3), int(size/2-2):int(size/2+3)])

    # sub_gal[[0,4,12,20,24]] = 0
    # plt.imshow(sub_gal.reshape((5,5)))
    # plt.show()
    # plt.close()

    z = numpy.delete(sub_gal, [0, 4, 12, 20, 24])

    a1, a2 = curve_fit(tool_box.fxy, x_d, z)
    snr = numpy.sqrt(10**a1[0]/numpy.sum(gal_p*rim)*num)
    flux_.append(snr)
flux_ = numpy.array(flux_)
plt.scatter(flux_, flux)
plt.show()
plt.hist(flux_-flux, 20)
plt.show()

plt.imshow(gals[10])
plt.show()