import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
# path.append('%s/work/mylib/' % my_home)
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import h5py
import numpy
import matplotlib.pyplot as plt
import tool_box


zbin_num = 6#int(argv[1])
theta_bin_num = 5#int(argv[2])
resample_num = 200#int(argv[3])

# discard the first bin
discard_bins = [0]

pts_num = int(theta_bin_num * (zbin_num ** 2 + zbin_num) / 2)
data_path = "E:/works/correlation/CFHT/cut_2.5/smooth"
pic_nm_p = data_path + "/xi_plus_mcmc_result_%d_compare.png" % resample_num
pic_nm_m = data_path + "/xi_minus_mcmc_result_%d_compare.png" % resample_num
pk_line_label = "MCMC:\n$As = 1.0\\times10^{-9}$\n$\Omega_m=0.27$"

pk_lines_tag = 0
if os.path.exists(data_path + "/mcmc_best_pk.npz"):
    npz = numpy.load(data_path + "/mcmc_best_pk.npz")
    xi_p_theoretical_lines = npz["arr_0"]
    xi_m_theoretical_lines = npz["arr_1"]
    xi_theta = npz["arr_2"]

    pk_lines_tag = 1
    print("Find Pk lines")

expo_type = ["diff_expo","same_expo"]

datas = []
cov = []
for ii in range(2):

    # pic_nm_p = data_path + "/chi_plus_result_%d_%s.png"%(resample_num,expo_type[ii])
    # pic_nm_m = data_path + "/chi_minus_result_%d_%s.png"%(resample_num,expo_type[ii])
    result_npz = data_path + "/result_cache_%d_%s.npz"%(resample_num,expo_type[ii])
    file_path = data_path + "/result_%d_%s.hdf5"%(resample_num,expo_type[ii])

    # print(list(h5f["/0"].keys()))

    results = tool_box.get_result_data(file_path, theta_bin_num,zbin_num, resample_num)
    theta, xi_p, xi_p_sig, xi_p_sub, xi_m, xi_m_sig, xi_m_sub = results

    use_cols = tool_box.get_zbin_mask(zbin_num, theta_bin_num, discard_bins)
    idx = use_cols > 0
    cov_p, inv_cov_p, cov_m, inv_cov_m, cov_pm, inv_cov_pm = tool_box.get_cov(xi_p_sub[idx], xi_m_sub[idx])
    n1, n2 = xi_p[idx].shape[0], xi_m[idx].shape[0]
    xi_pm = numpy.zeros((n1+n2, ))
    xi_pm[:n1] = xi_p[idx]
    xi_pm[n1:] = xi_m[idx]

    datas.append([theta, xi_p, xi_p_sig, xi_m, xi_m_sig, use_cols])

    cov.append([[cov_p, inv_cov_p], [cov_m, inv_cov_m], [cov_pm, inv_cov_pm]])

    numpy.savez(result_npz, theta[idx],
                xi_p[idx], xi_p_sig[idx], cov_p, inv_cov_p,
                xi_m[idx], xi_m_sig[idx], cov_m, inv_cov_m,
                xi_pm, cov_pm, inv_cov_pm)
    print(result_npz)

for i in range(len(cov)):

    img = Image_Plot(fig_x=4, fig_y=3, xpad=0.15, ypad=0.15, axis_linewidth=2.5,
                     plt_line_width=3, legend_size=25, xy_tick_size=25)
    img.subplots(3, 3)
    img.set_style()

    for j in range(3):
        for k in range(3):
            if k == 2:
                arr = cov[i][j][0]
                cov_scale = numpy.zeros_like(arr)
                y, x = arr.shape
                print(y,x)
                for p in range(y):
                    for q in range(x):
                        cov_scale[p, q] = arr[p, q] / numpy.sqrt(arr[p, p] * arr[q, q])

                fig = img.axs[j][k].imshow(cov_scale)
            else:
                fig = img.axs[j][k].imshow(cov[i][j][k])

            img.figure.colorbar(fig, ax=img.axs[j][k])
            img.del_axis(j,k,[0,1])

    img.axs[0][0].set_title("%s $\\xi_{+}$"%expo_type[i])
    img.axs[0][1].set_title("%s inv $\\xi_{+}$"%expo_type[i])
    img.axs[0][2].set_title("%s normalized $\\xi_{+}$"%expo_type[i])
    img.axs[1][0].set_title("%s $\\xi_{-}$" % expo_type[i])
    img.axs[1][1].set_title("%s inv $\\xi_{-}$" % expo_type[i])
    img.axs[1][2].set_title("%s normalized $\\xi_{-}$" % expo_type[i])
    img.axs[2][0].set_title("%s $\\xi_{+}+\\xi_{-}$" % expo_type[i])
    img.axs[2][1].set_title("%s inv $\\xi_{+}+\\xi_{-}$" % expo_type[i])
    img.axs[2][2].set_title("%s normalized $\\xi_{+}+\\xi_{-}$" % expo_type[i])

    img.save_img(data_path + "/cov_matrix_%s.png"%expo_type[i])
    img.close_img()


# plot xi_plus
img = Image_Plot(fig_x=4, fig_y=3, xpad=0, ypad=0, axis_linewidth=2.5, plt_line_width=3, legend_size=35,
                 xy_tick_size=25)
img.subplots(zbin_num, zbin_num)
img.set_style()


for ii in range(2):

    theta, xi_p, xi_p_sig = datas[ii][:3]
    use_cols = datas[ii][-1]

    if ii == 0:
        for i in range(zbin_num):
            for j in range(zbin_num-i,zbin_num):
                img.figure.delaxes(img.axs[i][j])

    img.axis_type(0,"major",10)
    img.axis_type(1,"major",10)
    img.axis_type(0,"minor",5)
    img.axis_type(1,"minor",5)

    tag = 0
    legend_tag = 0
    for i in range(zbin_num):
        for j in range(zbin_num):
            img_row = zbin_num - j - 1
            img_col = i
            st, ed = int(tag * theta_bin_num), int((tag + 1) * theta_bin_num)

            if j >= i:
                # print(theta[st:ed])
                alpha = 0.4
                if use_cols[st:ed].sum() > 1:
                    alpha = 1

                img.axs[img_row][img_col].errorbar(theta[st:ed]+theta[st:ed]*0.1*ii, xi_p[st:ed],yerr=xi_p_sig[st:ed],
                                                   lw=img.plt_line_width, elinewidth=img.plt_line_width, c="C%d"%ii,
                                                   marker="o",fmt=" ",mfc="none",ms=12,
                                                   label="$\\xi_{+}$ %s"%expo_type[ii],capsize=5, alpha=alpha)
                if use_cols[st:ed].sum() > 1 and legend_tag == 0:
                    legend_tag = 1
                    img.axs[img_row][img_col].legend(loc="lower left", bbox_to_anchor=(4, 1),
                                                     fancybox=False, fontsize=img.legend_size)

                img.axs_text(img_row, img_col, 0.83, 0.7, "%d-%d" % (i + 1, j + 1),
                             text_fontsize=img.legend_size-3, text_color="k")

                if ii == 0 and pk_lines_tag == 1:

                    img.axs[img_row][img_col].plot(xi_theta[tag], xi_p_theoretical_lines[tag], c="k",ls="dashdot", label=pk_line_label)

                img.axs[img_row][img_col].set_yscale("log")
                img.axs[img_row][img_col].set_xscale("log")

                img.axs[img_row][img_col].set_xlim(0.8, 60)
                img.set_ticklabel_str(img_row, img_col, 1,[1,5,10,20,40], ["$1$","$5$","$10$","$20$","$40$"])

                img.axs[img_row][img_col].set_ylim(4 * 10 ** (-7), 7 * 10 ** (-4))

                if img_col != 0:
                    img.axs[img_row][img_col].set_yticklabels([])

                tag += 1

img.save_img(pic_nm_p)
img.close_img()
print(pic_nm_p)
# img.show_img()


# plot xi_minus
img = Image_Plot(fig_x=4, fig_y=3, xpad=0, ypad=0, axis_linewidth=2.5, plt_line_width=3, legend_size=35,
                 xy_tick_size=25)
img.subplots(zbin_num, zbin_num)
img.set_style()
pic_nm = data_path + "/xi_minus_result_%d_compare.png" % resample_num

for ii in range(2):

    theta, xi_m, xi_m_sig = datas[ii][0],datas[ii][3],datas[ii][4]
    use_cols = datas[ii][-1]
    if ii == 0:
        for i in range(zbin_num):
            for j in range(zbin_num-i,zbin_num):
                img.figure.delaxes(img.axs[i][j])

    img.axis_type(0,"major",10)
    img.axis_type(1,"major",10)
    img.axis_type(0,"minor",5)
    img.axis_type(1,"minor",5)

    tag = 0
    legend_tag = 0
    for i in range(zbin_num):
        for j in range(zbin_num):
            img_row = zbin_num - j - 1
            img_col = i
            st, ed = int(tag * theta_bin_num), int((tag + 1) * theta_bin_num)
            if j >= i:
                # print(theta[st:ed])
                alpha = 0.4
                if use_cols[st:ed].sum() > 1:
                    alpha = 1
                img.axs[img_row][img_col].errorbar(theta[st:ed]+theta[st:ed]*0.1*ii, xi_m[st:ed],yerr=xi_m_sig[st:ed],
                                                   lw=img.plt_line_width, elinewidth=img.plt_line_width, c="C%d"%ii,
                                                   marker="o",fmt=" ",mfc="none",ms=12,
                                                   label="$\\xi_{+}$ %s"%expo_type[ii],capsize=5,alpha=alpha)
                if use_cols[st:ed].sum() > 1 and legend_tag == 0:
                    legend_tag = 1
                    img.axs[img_row][img_col].legend(loc="lower left", bbox_to_anchor=(4, 1),
                                                     fancybox=False, fontsize=img.legend_size)

                img.axs_text(img_row, img_col, 0.83, 0.7, "%d-%d" % (i + 1, j + 1),
                             text_fontsize=img.legend_size-3, text_color="k")

                if ii == 0 and pk_lines_tag == 1:
                    img.axs[img_row][img_col].plot(xi_theta[tag], xi_m_theoretical_lines[tag], c="k",ls="dashdot", label=pk_line_label)

                img.axs[img_row][img_col].set_yscale("log")
                img.axs[img_row][img_col].set_xscale("log")

                img.axs[img_row][img_col].set_xlim(0.8, 60)
                img.set_ticklabel_str(img_row, img_col, 1,[1,5,10,20,40], ["$1$","$5$","$10$","$20$","$40$"])

                img.axs[img_row][img_col].set_ylim(4 * 10 ** (-8), 7 * 10 ** (-5))

                if img_col != 0:
                    img.axs[img_row][img_col].set_yticklabels([])

                tag += 1

img.save_img(pic_nm_m)
img.close_img()
print(pic_nm_m)


