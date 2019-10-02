import numpy
from sys import path
path.append("D:/GitHub/astrophy-research/mylib/")
path.append("D:/GitHub/astrophy-research/multi_shear_detect")
import tool_box
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import gauss_fit_fun
import h5py


total_path = "E:/works/multi_shear/multi_shear_fit/4/"

for i in range(40):

    h5f = h5py.File(total_path + "cache_%d.hdf5"%i, "r")

    err_1s = h5f["chisq_1"].value
    para_1s = h5f["/para_1"].value

    err_2s = h5f["chisq_2"].value
    para_2s = h5f["/para_2"].value

    signal = h5f["/signal"].value
    sigma = h5f["/sigma"].value
    x_fit = h5f["/x_fit"].value
    y_fit = h5f["/y_fit"].value

    h5f.close()

    err_1s.shape = (len(err_1s), 1)
    err_2s.shape = (len(err_2s), 1)
    # print(para_1s.shape,para_2s.shape, len(para_1s))
    # para_1s.shape = (len(para_1s), para_1s.shape[1])
    # para_2s.shape = (len(para_2s), para_2s.shape[1])

    if i == 0:
        chisq_1 = err_1s
        chisq_2 = err_2s
        all_para_1 = para_1s
        all_para_2 = para_2s
    else:
        chisq_1 = numpy.row_stack((chisq_1, err_1s))
        chisq_2 = numpy.row_stack((chisq_2, err_2s))
        all_para_1 = numpy.row_stack((all_para_1, para_1s))
        all_para_2 = numpy.row_stack((all_para_2, para_2s))
    continue

    print("\n%d"%i)
    print("%d, %.4f, %.4f"%(err_1s.shape[0], err_1s.min(), err_1s.max()))
    print("%d, %.4f, %.4f"%(err_2s.shape[0], err_2s.min(), err_2s.max()))

    idx = err_1s <= 10
    if idx.sum() > 0:
        print(idx.sum(), err_1s[idx].min(), err_1s[idx].max())
    else:
        print(0)
    idx = err_2s <= 10
    if idx.sum() > 0:
        print(idx.sum(), err_2s[idx].min(), err_2s[idx].max())
    else:
        print(0)
    #
    # idx = err_1s > 0.34
    # print(idx.sum(), err_1s[idx].min())
    # idx = err_2s > 0.34
    # print(idx.sum(), err_2s[idx].min())
    #
    #
    # histogram of \chi^2
    img = Image_Plot()
    img.subplots(1,2)
    img.axs[0][0].hist(err_1s,100)
    img.axs[0][1].hist(err_2s,100)
    img.set_label(0,0,1,"$\chi^2$")
    img.set_label(0,1,1,"$\chi^2$")
    img.set_label(0,0,0,"Num")
    img.set_label(0,1,0,"Num")
    img.axs[0][0].set_yscale("log")
    img.axs[0][1].set_yscale("log")
    img.axs[0][0].set_xscale("log")
    img.axs[0][1].set_xscale("log")
    img.save_img(total_path + "chisq_%d.png"%i)
    img.close_img()
    # img.show_img()
    #
    #
    # plot the good fitting result
    # the input one
    fit_true = gauss_fit_fun.gauss_2_coeff(x_fit, 0.5, signal[0], sigma[0], 0.5, signal[1], sigma[0])
    fit_true1 = gauss_fit_fun.gauss_coeff(x_fit, 0.5, signal[0], sigma[0])
    fit_true2 = gauss_fit_fun.gauss_coeff(x_fit, 0.5, signal[1], sigma[0])
    chisq_true = gauss_fit_fun.fit_chisq_2_coeff([0.5, signal[0], sigma[0], 0.5, signal[1], sigma[0]], x_fit, y_fit)
    print("True chi square: ", chisq_true)

    # one gaussian's fitting
    idx = err_1s == err_1s.min()
    para1 = para_1s[idx][0]
    print("The minimum: %.5f"%err_1s[idx][0], para1)
    fit_s = gauss_fit_fun.gauss_coeff(x_fit, para1[0], para1[1], para1[2])

    # two gaussian's fitting
    idx = err_2s == err_2s.min()
    num1, mu1, sig1, num2, mu2, sig2 = para_2s[idx][0]
    if idx.sum() > 1:
        print("The minimum: %.5f"%err_2s[idx][0], num1, mu1, sig1, num2, mu2, sig2)
    else:
        print("The minimum: %.5f"%err_2s[idx], num1, mu1, sig1, num2, mu2, sig2)
    fit_d = gauss_fit_fun.gauss_2_coeff(x_fit, num1, mu1, sig1, num2, mu2, sig2)
    fit_d1 = gauss_fit_fun.gauss_coeff(x_fit, num1, mu1, sig1)
    fit_d2 = gauss_fit_fun.gauss_coeff(x_fit, num2, mu2, sig2)


    img = Image_Plot(plt_line_width=3, fig_x=12, fig_y=9)
    img.subplots(1, 2)

    for j in range(2):
        img.axs[0][j].plot(x_fit, fit_true, linewidth=img.plt_line_width+1, linestyle="-",c="C0",
                           label="True f1 + f2, $\chi^2=%.3f$"%chisq_true)
        img.axs[0][j].plot(x_fit, fit_true1, linewidth=img.plt_line_width+1, linestyle="dotted",c="C5",label="True f1")
        img.axs[0][j].plot(x_fit, fit_true2, linewidth=img.plt_line_width+1, linestyle="dotted",c="C6",label="True f2")
        img.axs[0][j].scatter(x_fit, y_fit, c="black", s=20, label="data")


    img.axs[0][0].plot(x_fit, fit_s, linewidth=img.plt_line_width, c="C1", linestyle="--",
                       label="$\chi^2=%.3f$ \n%.3f*N($\mu=%.4f$, $\sigma=%.4f$)"
                             % (err_1s.min(), para1[0], para1[1], para1[2]))

    img.axs[0][1].plot(x_fit, fit_d, linewidth=img.plt_line_width, c="C1", linestyle="--",
                       label="f1 + f2, $\chi^2=%.3f$" % err_2s.min())
    img.axs[0][1].plot(x_fit, fit_d1, linewidth=img.plt_line_width, c="C2",linestyle="-.",
                       label="f1 %.3f*N($\mu=%.4f$, $\sigma=%.4f$)" % (num1, mu1, sig1))
    img.axs[0][1].plot(x_fit, fit_d2, linewidth=img.plt_line_width, c="C4",linestyle="-.",
                       label="f2 %.3f*N($\mu=%.4f$, $\sigma=%.4f$)" % (num2, mu2, sig2))
    for j in range(2):
        img.set_label(0,j,0,"P(G)")
        img.set_label(0,j,1,"G")
        img.axs[0][j].legend(fontsize=img.legend_size,loc="upper left",frameon=False)
    img.save_img(total_path + "fit_%d.png"%i)
    img.close_img()
    # img.show_img()

# numpy.savez(total_path + "chisq_cache.npz", chisq_1, chisq_2)

chisq_1.shape = (chisq_1.shape[0],)
chisq_2.shape = (chisq_2.shape[0],)
print(chisq_1.shape, all_para_1.shape, all_para_2.shape)

idx1 = chisq_2 <= 0.02
idx21 = chisq_2 >= 690
idx22 = chisq_2 <= 700
idx2 = idx21 & idx22
idx3 = chisq_2 >= 700

# img = Image_Plot(fig_x=9, fig_y=6)
# img.subplots(1,3)
# img.axs[0][0].hist(chisq_2[idx1], 100)
# img.axs[0][1].hist(chisq_2[idx2], 100)
# img.axs[0][2].hist(chisq_2[idx3], 100)
# # img.axs[0][2].hist(chisq_1, 100)
# img.axs[0][0].set_yscale("log")
# img.axs[0][1].set_yscale("log")
# img.axs[0][2].set_yscale("log")
# img.show_img()
#
titles = ["weight_1", "mu_1", "sgima_1","weight_2", "mu_2", "sgima_2"]
# img = Image_Plot()
# img.subplots(3,6)
# index = [idx1, idx2, idx3]
# for i in range(3):
#     for j in range(6):
#         plt_data = all_para_2[:, j][index[i]]
#         img.axs[i][j].hist(plt_data, 10)
#         img.axs[i][j].set_title("%.4f <= %s <=%.4f"%(plt_data.min(), titles[j],plt_data.max()))
#         img.axs[i][j].set_yscale("log")
#
# img.save_img(total_path + "para_hist.png")
# img.show_img()

img = Image_Plot(fig_x=6,fig_y=4)
img.subplots(5,5)
for i in range(5):
    data_1 = all_para_2[:, i][idx1]
    for j in range(i+1, 6):
        print(titles[i], titles[j])
        data_2 = all_para_2[:, j][idx1]
        chisq = chisq_2[idx]
        norm = img.figure.Normalize(vmin=numpy.min(chisq), vmax=numpy.max(chisq))
        cmap = img.figure.get_cmap('YlOrRd')
        img.axs[i][j-i-1].scatter(data_2, data_1, s=10)
        img.set_label(i,j-i-1, 0, titles[i])
        img.set_label(i,j-i-1, 1, titles[j])

for i in range(5):
    for j in range(4,4-i,-1):
        print(i,j)
        img.del_tick(i,j,[0,1], [0,1,2,3])
img.subimg_adjust(0.27,0.27)
img.save_img(total_path + "para_contour_.png")
img.show_img()
# print(idx1.sum(), idx2.sum())
# print(chisq_1[idx1].min(), chisq_1[idx1].max())
# print(chisq_1[idx2].min(), chisq_1[idx2].max())




# print(chisq_2.min())
# idx1 = chisq_2 <= 0.03
# idx21 = chisq_2 > 0.03
# idx22 = chisq_2 < 0.2
# idx2 = idx21 & idx22
# idx3 = chisq_2 > 10
# img = Image_Plot()
# img.subplots(1, 1)
# img.axs[0][0].hist(chisq_2[idx1], 2)
# img.axs[0][0].hist(chisq_2[idx2], 10)
# img.axs[0][0].set_yscale("log")
# img.show_img()
#
# print(chisq_2.shape[0])
# print(chisq_2[idx1].min(), chisq_2[idx1].max())
# print(chisq_2[idx2].min(), chisq_2[idx2].max())
# print(chisq_2[idx3].min(), chisq_2[idx3].max())
# print(idx1.sum(), idx2.sum(), idx3.sum())

#
#
# if "7" in total_path:
#     print(chisq_1.shape)
#     idx1 = chisq_1 <= 10
#     idx21 = chisq_1 > 350
#     idx22 = chisq_1 < 370
#     idx2 = idx21 & idx22
#     idx31 = chisq_1 > 1000
#
#     print(chisq_1.shape[0])
#     print(idx2.sum(), idx31.sum())
#     # print(chisq_1[idx1].min(),chisq_1[idx1].max())
#     print(chisq_1[idx2].min(), chisq_1[idx2].max())
#     print(chisq_1[idx31].min(), chisq_1[idx31].max())
#
#
#
#     idx1 = chisq_2 <= 10
#     idx21 = chisq_2 > 350
#     idx22 = chisq_2 < 380
#     idx2 = idx21 & idx22
#     idx31 = chisq_2 > 1000
#     # print(chisq_2)
#     print(chisq_2.shape[0])
#     print(idx1.sum(), idx2.sum(), idx31.sum())
#     print(chisq_2[idx1].min(),chisq_2[idx1].max())
#     print(chisq_2[idx2].min(), chisq_2[idx2].max())
#     print(chisq_2[idx31].min(), chisq_2[idx31].max())
#
