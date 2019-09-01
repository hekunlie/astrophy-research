import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
import h5py
from Fourier_Quad import Fourier_Quad

fq = Fourier_Quad(4,14)
tag = [0,1,2,5]

# for i in tag:
#     h5f = h5py.File("D:/result_%d.hdf5"%i,"r")
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


h5f = h5py.File("D:/result_0.hdf5","r")
data = h5f["/data"].value
mg1 = data[:,0]
mg2 = data[:,1]
mnu1 = data[:,2] + data[:,3]
mnu2 = data[:,2] - data[:,3]
img = Image_Plot()
img.subplots(1,1)
gh = numpy.linspace(-0.05, 0.05, 3)

a1, a2 = -5.e6, 5.e6
idx1 = mg1 >= a1
idx2 = mg1 <= a2

bin_num = 200

nums, bins = numpy.histogram(mg1[idx1&idx2], 200)[:2]
bins = bins[:-1]
npw = numpy.where(nums == nums.max())[0][0]
print(npw,nums.shape,bins.shape)

length = min(npw-0, bin_num-1-npw)

line_w = 3
# total PDF
img.axs[0][0].scatter(bins, nums, label="PDF",c="C9", s=15, marker="p")

perc = [1./10, 2./10, 3./10, 5./10, 7./10, 1]
for i in range(len(perc)):
    # fit with full range
    st, ed = npw - int(perc[i]*length), npw + int(perc[i]*length)
    fit_bin = bins[st:ed]
    fit_nums = nums[st:ed]

    para = tool_box.gauss_fit([fit_bin], fit_nums)
    a, b, c = para[1][0], para[0][0], para[0][1]

    fx = a*numpy.exp(-(fit_bin-b)**2/2/abs(c))
    img.axs[0][0].plot(fit_bin, fx, linewidth=line_w, linestyle="--")

    # # fit with central part
    # bin_st, bin_ed = 90, 110
    # fit_bin = bins[bin_st:bin_ed]
    # para = tool_box.gauss_fit([fit_bin], nums[bin_st:bin_ed])
    # a, b, c = para[1][0], para[0][0], para[0][1]
    # fx = a*numpy.exp(-(fit_bin-b)**2/2/abs(c))
    # img.axs[0][0].plot(fit_bin, fx, linewidth=line_w, linestyle="--",c="C3", label="Fit PDF 2")
    #
    # # fit with central part
    # bin_st, bin_ed = 50, 150
    # fit_bin = bins[bin_st:bin_ed]
    # para = tool_box.gauss_fit([fit_bin], nums[bin_st:bin_ed])
    # a, b, c = para[1][0], para[0][0], para[0][1]
    # fx = a*numpy.exp(-(fit_bin-b)**2/2/abs(c))
    # img.axs[0][0].plot(fit_bin, fx, linewidth=line_w, linestyle="--",c="C4", label="Fit PDF 3")

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
img.axs[0][0].set_xlim(a1,a2)
img.show_img()