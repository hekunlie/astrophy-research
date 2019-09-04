import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib/")
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
import h5py
from Fourier_Quad import Fourier_Quad
import scipy


def mix_gauss(p, x):
    a1,a2,a3, b1, b2, b3 = p
    f1 = a1/numpy.sqrt(2*numpy.pi)/a2*numpy.exp(-(x-a3)**2/2/a2**2)
    f2 = b1/numpy.sqrt(2*numpy.pi)/b2*numpy.exp(-(x-b3)**2/2/b2**2)
    return f1 + f2

def error(p, x, y):
    return mix_gauss(p, x) - y


fq = Fourier_Quad(4,14)
tag = [0,1,2,3]

# for i in tag:
#     h5f = h5py.File("F:/works/multi_shear_detect/sig_50/result_%d.hdf5"%i,"r")
#     data = h5f["/data"].value
#     h5f.close()
#
#     mg1 = data[:,0]
#     mg2 = data[:,1]
#     mnu1 = data[:,2] + data[:,3]
#     mnu2 = data[:,2] - data[:,3]
#
#     gh, gh_sig = fq.find_shear(mg1, mnu1, 8)[:2]
#
#     print(gh, gh_sig)


h5f = h5py.File("F:/works/multi_shear_detect/sig_10/result_0.hdf5","r")
data = h5f["/data"].value
mg1 = data[:,0]
mg2 = data[:,1]
mnu1 = data[:,2] + data[:,3]
mnu2 = data[:,2] - data[:,3]



gh = numpy.linspace(-0.05, 0.05, 3)

a1, a2 = -9.e6, 9.e6
idx1 = mg1 >= a1
idx2 = mg1 <= a2

bin_num = 500

xscale = 1./10000
nums, bins = numpy.histogram(mg1[idx1&idx2], bin_num)[:2]
nums = nums/nums.sum()
bins = bins[:-1]*xscale
npw = numpy.where(nums == nums.max())[0][0]
bin_width = bins[2] - bins[1]
print(npw, nums.shape, bins.shape)

length = min(npw-0, bin_num-1-npw)

line_w = 3


img = Image_Plot()
img.subplots(1,2)
# total PDF
# img.axs[0][0].bar(bins, nums, width=bin_width, alpha=0.5,label="PDF",color="none",edgecolor="C9")
img.axs[0][0].fill_between(bins, 0, nums, alpha=0.5,label="PDF",color="C9")
perc = [1./10, 2./10, 3./10, 5./10, 7./10, 1]


coeff_0 = [10, 10, 0, 10, 10, 0]
res = scipy.optimize.least_squares(error, coeff_0, args=(bins, nums), method="trf")
fx_fit = mix_gauss(res.x,bins)
# img.axs[0][1].scatter(bins, nums, label="PDF",c="C9", s=15, marker="p")
img.axs[0][1].fill_between(bins, 0, nums, alpha=0.5, label="PDF",color="C9")
img.axs[0][1].plot(bins, fx_fit, linewidth=line_w, linestyle="--", c="C1", label="Fitted by two Gaussian's")

for i in range(len(perc)):
    # fit with full range
    st, ed = npw - int(perc[i]*length), npw + int(perc[i]*length)
    fit_bin = bins[st:ed]
    fit_nums = nums[st:ed]

    para = tool_box.gauss_fit([fit_bin], fit_nums)
    a, b, c = para[1][0], para[0][0], para[0][1]

    fx = a*numpy.exp(-(fit_bin-b)**2/2/abs(c))
    img.axs[0][0].plot(fit_bin, fx, linewidth=line_w, linestyle="--", c="C%d"%i, label="%.2f"%perc[i])


# bin_num = 8
# bins = fq.set_bin(mg1, bin_num, 100)
#
# bin_num2 = int(bin_num * 0.5)
# inverse = range(int(bin_num / 2 - 1), -1, -1)
# for ig in gh:
#     G_hat = mg1 - ig*mnu1
#     chisq = fq.get_chisq(mg1, mnu1, ig, bins, bin_num2, inverse, 0)
#     img.axs[0][0].hist(G_hat, 50000, histtype="step", label="$\hat g$=%.3f,$\chi^2$=%.2f"%(ig,chisq))
img.axs[0][0].legend()
img.axs[0][1].legend()
img.axs[0][0].set_xlim(a1*xscale,a2*xscale)
img.show_img()