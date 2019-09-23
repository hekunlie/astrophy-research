import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import h5py
from scipy import optimize



def gauss(x, mu, sig):
    return numpy.exp(-(x-mu)**2/2/sig/sig)

def fit_pdf(x, w1, mu1, sig1, w2, mu2, sig2):
    return w1*numpy.exp(-(x-mu1)**3/2/sig1) + w2*x*numpy.exp(-(x-mu2)**2/2/sig2)


h5f = h5py.File("E:/CFHT_W1.hdf5","r")
data = h5f["/data"].value
z = data[:,4]
zmin = data[:,5]
zmax = data[:,6]
zodd = data[:,7]
print(data[0])

img = Image_Plot()
img.subplots(2,2)
num, bins = img.axs[0][0].hist(z, 100, density=True)[:2]
x_fit = (bins[1:] + bins[:bins.shape[0]-1])/2
y_fit = num
res = optimize.curve_fit(fit_pdf,x_fit,y_fit, bounds=[0,[10,2,4,10,5,4]])[0]
w1, mu1, sig1, w2, mu2, sig2 = res
img.axs[0][0].plot(x_fit, fit_pdf(x_fit,w1, mu1, sig1, w2, mu2, sig2))
print(res)
img.axs[0][1].hist((zmax-zmin)/2, 50, density=True)
img.axs[1][0].hist(zodd, 50, density=True)
img.axs[1][1].scatter(z[:200000], zodd[:200000], s=3)
img.show_img()