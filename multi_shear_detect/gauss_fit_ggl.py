from sys import path
path.append("/home/hklee/work/mylib/")
path.append("D:/GitHub/astrophy-research/mylib/")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import h5py
import numpy
import gauss_fit_fun
import scipy
import tool_box


coeff = 554.682135528

fq = Fourier_Quad(12,123)

total_path = "E:/GGL_FQ/"
# signal = numpy.zeros((2,10))

npz_cache = numpy.load("%scache.npz"%total_path)
signal = npz_cache["arr_0"]
print(signal.shape)
r = tool_box.set_bin_log(0.04,15,13)[:10]
# img = Image_Plot(plt_line_width=2)
# img.subplots(1, 1)
# img.set_style()
# img.axs[0][0].errorbar(r,signal[0],signal[1],capsize=3)
# img.axs[0][0].set_xscale("log")
# img.axs[0][0].set_yscale("log")
# img.set_label(0,0,0,"$\Delta\Sigma \; [\\rm{h \cdot M_{\odot}} \cdot \\rm{pc^{-2}}]$")
# img.set_label(0,0,1,"$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$")
# img.save_img(total_path + "signal.png")
# img.show_img()
# exit()
nums = []
for ig in range(10):
    h5f = h5py.File("%sradius_%d.hdf5"%(total_path,ig),"r")
    data_col = h5f["/data_col"].value
    pair_count = h5f["/pair_count"].value
    data_ori = h5f["/pair_data_0"].value
    radius_bin = h5f["/radius_bin"].value
    h5f.close()
    nums.append(data_ori.shape[0])

    redshift = data_ori[:, 6]
    idx1 = redshift >= 0
    idx2 = redshift <= 100

    idx = idx1 & idx2
    data = data_ori[idx]

    sigma = data[:,4]*coeff
    mg1 = data[:,0]
    mg2 = data[:,1]
    mnu1 = data[:,2]
    mnu2 = data[:,3]
    radius_sep = data[:, 5]
    redshift = data[:, 6]
    del_sigma, sigma_sig = fq.find_shear(mg1*sigma, mnu1, 8, left=0, right=120)[:2]
    # del_sigma, sigma_sig = signal[0,ig],signal[1,ig]
    print(del_sigma, sigma_sig)

    bin_num = 11
    gh = numpy.linspace(-10,10, bin_num)

    x_fit = (gh[:bin_num - 1] + gh[1:]) / 2
    dg = gh[1] - gh[0]
    num_move = gauss_fit_fun.get_flow(mg1*sigma, mnu1, gh)
    y_fit = num_move/num_move.sum()/dg

    # fit single gauss
    # fit_res = scipy.optimize.curve_fit(gauss_fit_fun.gauss_coeff, x_fit, y_fit*1000, method="trf")[0]
    # fx_fit = gauss_fit_fun.gauss_coeff(x_fit, fit_res[0]/1000,fit_res[1],fit_res[2])
    # fx_fit = gauss_fit_fun.gauss_coeff(x_fit, 1, 0, fit_res[2])

    res = tool_box.gauss_fit([x_fit],y_fit)
    fit_res = [res[1][0],res[0][0],numpy.sqrt(res[0][1])]
    fx_fit = fit_res[0]*numpy.exp(-(x_fit - fit_res[1])**2/2/fit_res[2]**2)
    print(res)
    fx_ratio = y_fit / fx_fit - 1

    img = Image_Plot(plt_line_width=2)
    img.subplots(1, 3)
    img.set_style()

    img.axs[0][0].scatter(x_fit, y_fit, label="Num: %d,\n$\Delta \Sigma=%.3f(%.3f)$"%(len(redshift),del_sigma, sigma_sig))
    img.axs[0][0].plot(x_fit, fx_fit, c="C1",label="$%.5fN(%.3f,%.3f)$"%(fit_res[0],fit_res[1],fit_res[2]),linewidth=img.plt_line_width)

    dy = (y_fit.max() - y_fit.min())*0.1
    img.axs[0][0].set_ylim(y_fit.min()-dy, y_fit.max()+dy)
    img.set_label(0,0,0,"$P(\hat{\Delta\Sigma})$")
    img.set_label(0,0,1,"$\hat{\Delta\Sigma}$")
    img.axs[0][0].legend(fontsize=img.legend_size,loc="lower center")


    img.axs[0][1].plot(x_fit, fx_ratio,)
    img.set_label(0,1,0,"$\\frac{P(\Delta\hat\Sigma)}{P(\Delta\hat \Sigma)_{fit}} - 1$")
    img.set_label(0,1,1,"$\hat{\Delta\Sigma}$")
    # img.axs[0][1].legend(fontsize=img.legend_size)

    if ig <= 5:
        hist_bin_num = 501
    else:
        hist_bin_num = 2001

    hist_bins = numpy.linspace(-1.e6,1.e6,hist_bin_num)
    img.axs[0][2].hist(mg1*sigma,hist_bins,histtype="step",label="$\hat{\Delta\Sigma}=0$",linewidth=img.plt_line_width)
    img.axs[0][2].hist(mg1*sigma-mnu1*gh[0],hist_bins,histtype="step",label="$\hat{\Delta\Sigma}=%.2f$"%gh[0],linewidth=img.plt_line_width)
    img.axs[0][2].hist(mg1*sigma-mnu1*gh[-1],hist_bins,histtype="step",label="$\hat{\Delta\Sigma}=%.2f$"%gh[-1],linewidth=img.plt_line_width)

    img.set_label(0,2,0,"$Num$")
    img.set_label(0,2,1,"$G_1\\frac{D_S}{D_S D_{LS}} - \hat{\Delta\Sigma}(N+U)$")

    img.axs[0][2].legend(fontsize=img.legend_size)
    img.axs[0][2].set_xlim(-1.e5,1.e5)
    img.subimg_adjust(0,0.27)
    img.save_img(total_path + "test/%d_1.png"%ig)
    # img.show_img()
    img.close_img()

img = Image_Plot(plt_line_width=2)
img.subplots(1, 1)
img.set_style()
img.axs[0][0].scatter(r, nums)
img.axs[0][0].set_xscale("log")
img.axs[0][0].set_yscale("log")
img.set_label(0,0,0,"Pair Num")
img.set_label(0,0,1,"$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$")
img.save_img(total_path + "num.png")
img.show_img()