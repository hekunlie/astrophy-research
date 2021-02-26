from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box
import time
import numpy


chi_guess_num = 15
chi_guess = tool_box.set_bin_log(1./10**9,1./10**4,chi_guess_num)
print(chi_guess)

chi_measures = [numpy.zeros((2, chi_guess_num)),
                numpy.zeros((2, chi_guess_num))]

data_path = "D:/TEMP/chisq/2/narrow/even_pts"
data_type = ["noise_free", "noisy_cpp"]
titles = ["noise_free,<$g_1g_1^{\prime}$>, $1\\times 10^8$ gals",
          "noisy,<$g_1g_1^{\prime}$>, $1\\times 10^8$ gals"]

for j in range(2):
    img = Image_Plot(fig_x=14,fig_y=5,xpad=0.15,ypad=0.15)
    img.subplots(1,1)


    for i in range(10,7,-1):
    # for i in range(chi_guess_num-1,-1,-1):
        tag_1 = int(2*i)
        tag_2 = tag_1 + 1

        npz = numpy.load(data_path + "/chisq_%d_%d_%s_narrow.npz"%(tag_1, tag_2, data_type[j]))
        chi_guess_bin = npz["arr_0"]
        chisq_arr = npz["arr_1"][0]
        idx = chisq_arr <= 150
        ls = img.line_styles[i]

        coeff = tool_box.fit_1d(chi_guess_bin[idx], chisq_arr[idx], 2, "scipy")
        xi_err = numpy.sqrt(1 / 2. / coeff[2])
        xi = -coeff[1] / 2. / coeff[2]
        chi_measures[j][0,i] = xi
        chi_measures[j][1,i] = xi_err
        print(tag_1,tag_2,i,chi_guess[i],xi,xi_err)

        img.axs[0][0].plot(chi_guess_bin[idx], chisq_arr[idx], linewidth=2.2,
                           label="$\\xi^t=%.2e, \\xi^m=%.3e(%.2e)$"%(chi_guess[i],xi,xi_err), ls=ls[0],c=ls[1])
        ys = img.axs[0][0].set_ylim()
        img.axs[0][0].plot([chi_guess[i],chi_guess[i]],[ys[0],ys[1]],c=ls[1], label="True \\xi=%.3e$")
    img.axs[0][0].legend(ncol=2,fontsize=12,loc="lower center")
    ys = img.axs[0][0].set_ylim()
    img.axs[0][0].set_xscale("symlog")

    img.axs[0][0].set_ylim(1, 170)
    # img.axs[0][0].set_xlim(-4*10**(-6),8*10**(-6))
    img.axs[0][0].set_yscale("log")
    img.axs[0][0].set_title(titles[j])
    img.set_label(0, 0, 0, "$\\chi^2$")
    img.set_label(0, 0, 1, "$\\xi_1$ guess")
    img.set_ticklabel_str(0, 0, 1, [1.e-7,1.e-6, 1.e-5],
                          ["$10^{-%d}$" % i for i in range(7, 4, -1)])

    img.save_img(data_path + "/chisq_%s.png"%data_type[j])
    # img.show_img()
    img.close_img()
numpy.savez(data_path+"/xi_cache_%s.npz"%data_type[0],chi_measures[0])
numpy.savez(data_path+"/xi_cache_%s.npz"%data_type[1],chi_measures[1])

cs = ["C2","C1"]
img = Image_Plot(fig_x=5,fig_y=5,xpad=0.15,ypad=0.15)
img.subplots(1,1)
for i in range(1,-1,-1):
    img.axs[0][0].errorbar(chi_guess, chi_measures[i][0],chi_measures[i][1], marker="s",ms=3,c=cs[i], capsize=3, fmt=" ", label=titles[i])

img.axs[0][0].plot(chi_guess,chi_guess, c="k", ls="--", label="y=x")
img.axs[0][0].legend()
img.axs[0][0].set_xscale("log")
img.axs[0][0].set_yscale("log")
img.set_ticklabel_str(0,0,1,[1.e-9,1.e-8,1.e-7,1.e-6,1.e-5,1.e-4],["$10^{-%d}$"%i for i in range(9,3,-1)])
img.set_label(0,0,0,"true $\\xi_1$")
img.set_label(0,0,1,"measured $\\xi_1$")
img.axs[0][0].grid()
img.save_img(data_path + "/compare.png")
# img.show_img()