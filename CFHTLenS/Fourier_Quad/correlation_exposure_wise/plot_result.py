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

pts_num = int(theta_bin_num * (zbin_num ** 2 + zbin_num) / 2)
data_path = "E:/works/correlation/CFHT/cut_2.5/smooth"

pk_lines_tag = 0
if os.path.exists(data_path + "/tomo_pk.npz"):
    npz = numpy.load(data_path + "/tomo_pk.npz")
    pk_lines = npz["arr_6"]
    pk_theta = npz["arr_5"]

    pk_lines_tag = 1
    print("Find Pk lines")

# chi_plus
img = Image_Plot(fig_x=4, fig_y=3,xpad=0,ypad=0,axis_linewidth=2.5, plt_line_width=3, legend_size=35,xy_tick_size=25)
img.subplots(zbin_num, zbin_num)
img.set_style()
pic_nm_p = data_path + "/chi_plus_result_%d_compare.png" % resample_num

expo_type = ["diff_expo","same_expo"]

cov = []
for ii in range(2):

    # pic_nm_p = data_path + "/chi_plus_result_%d_%s.png"%(resample_num,expo_type[ii])
    # pic_nm_m = data_path + "/chi_minus_result_%d_%s.png"%(resample_num,expo_type[ii])
    result_npz = data_path + "/result_cache_%d_%s.npz"%(resample_num,expo_type[ii])
    file_path = data_path + "/result_%d_%s.hdf5"%(resample_num,expo_type[ii])

    # print(list(h5f["/0"].keys()))

    theta, xi_p, xi_p_sig, cov_p, xi_m, xi_m_sig, cov_m, xi_pm, cov_pm = tool_box.get_result_data(file_path, pts_num, resample_num)

    inv_cov_p = numpy.linalg.pinv(cov_p*(resample_num-1))
    inv_cov_m = numpy.linalg.pinv(cov_m*(resample_num-1))
    inv_cov_pm = numpy.linalg.pinv(cov_pm*(resample_num-1))

    cov.append([[cov_p, inv_cov_p], [cov_m, inv_cov_m], [cov_pm, inv_cov_pm]])

    numpy.savez(result_npz, theta,
                xi_p, xi_p_sig, cov_p, inv_cov_p,
                xi_m, xi_m_sig, cov_m, inv_cov_m,
                xi_pm, cov_pm, inv_cov_pm)

    if ii == 0:
        for i in range(zbin_num):
            for j in range(zbin_num-i,zbin_num):
                img.figure.delaxes(img.axs[i][j])

    img.axis_type(0,"major",10)
    img.axis_type(1,"major",10)
    img.axis_type(0,"minor",5)
    img.axis_type(1,"minor",5)

    tag = 0
    for i in range(zbin_num):
        for j in range(zbin_num):
            img_row = zbin_num - j - 1
            img_col = i
            if j >= i:
                st, ed = int(tag*theta_bin_num), int((tag+1)*theta_bin_num)
                # print(theta[st:ed])

                img.axs[img_row][img_col].errorbar(theta[0,st:ed]+theta[0,st:ed]*0.1*ii, xi_p[0,st:ed],yerr=xi_p_sig[0,st:ed],
                                                   lw=img.plt_line_width, elinewidth=img.plt_line_width, c="C%d"%ii,
                                                   marker="o",fmt=" ",mfc="none",ms=12,label="$\\xi_{+}$ %s"%expo_type[ii],capsize=5)
                if tag == 0:
                    img.axs[img_row][img_col].legend(loc="lower left",bbox_to_anchor=(4,1),
                                                     fancybox=False,fontsize=img.legend_size)

                img.axs_text(img_row, img_col, 0.8, 0.8, "%d-%d" % (i + 1, j + 1), text_fontsize=img.legend_size, text_color="k")

                if ii == 0 and pk_lines_tag == 1:
                    pk_line_label = "Planck2018:\n$\Omega_m=0.313$\n$\sigma_8=0.8097$"
                    img.axs[img_row][img_col].plot(pk_theta, pk_lines[tag], c="k",ls="dashdot", label=pk_line_label)


                # img.axs[img_row][img_col].set_xticks(fontsize=img.legend_size)
                # img.axs[img_row][img_col].set_yticks(fontsize=img.legend_size)
                img.axs[img_row][img_col].set_yscale("log")
                img.axs[img_row][img_col].set_xscale("log")

                img.axs[img_row][img_col].set_xlim(0.8, 60)
                img.set_ticklabel_str(img_row, img_col, 1,[1,5,10,20,40], ["$1$","$5$","$10$","$20$","$40$"])

                # img.axs[img_row][img_col].set_xlim(0.8, 200)
                # img.set_ticklabel_str(img_row, img_col, 1,[1,5,10,100], ["$1$","$5$","$10$","$100$"])
                img.axs[img_row][img_col].set_ylim(4 * 10 ** (-7), 7 * 10 ** (-4))

                # if tag not in [0, 6, 11, 15, 18, 20]:
                #     # img.del_axis(img_row,img_col,[1])
                #     img.axs[img_row][img_col].set_xticklabels([])
                if img_col != 0:
                    img.axs[img_row][img_col].set_yticklabels([])
                    # img.axs[img_row][img_col].set_yticks([])

                tag += 1



img.save_img(pic_nm_p)
img.close_img()
print(pic_nm_p)
# img.show_img()

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
                print(y, x)
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

exit()

# chi_minus
img = Image_Plot(fig_x=4, fig_y=3,xpad=0,ypad=0,axis_linewidth=2.5, plt_line_width=3, legend_size=25,xy_tick_size=25)
img.subplots(zbin_num, zbin_num)
img.set_style()

for i in range(zbin_num):
    for j in range(zbin_num-i,zbin_num):

        img.figure.delaxes(img.axs[i][j])

img.axis_type(0,"major",10)
img.axis_type(1,"major",10)
img.axis_type(0,"minor",5)
img.axis_type(1,"minor",5)

tag = 0
for i in range(zbin_num):
    for j in range(zbin_num):
        img_row = zbin_num - j - 1
        img_col = i
        if j >= i:
            st, ed = int(tag*theta_bin_num), int((tag+1)*theta_bin_num)
            # print(theta[st:ed])
            img.axs[img_row][img_col].errorbar(theta[0,st:ed], xi_m[0,st:ed],yerr=xi_m_sig[0,st:ed],
                                               lw=img.plt_line_width,elinewidth=img.plt_line_width,
                                               marker="o",fmt=" ",color="C1",mfc="none",ms=12,label="$\\xi_{-}$",capsize=5)
            if tag == 0:
                img.axs[img_row][img_col].legend(loc="lower left",bbox_to_anchor=(1.2,0.25),
                                                 fancybox=False,fontsize=img.legend_size)

            img.axs_text(img_row, img_col, 0.8, 0.8, "%d-%d" % (i + 1, j + 1), text_fontsize=img.legend_size, text_color="k")

            # img.axs[img_row][img_col].set_xticks(fontsize=img.legend_size)
            # img.axs[img_row][img_col].set_yticks(fontsize=img.legend_size)
            img.axs[img_row][img_col].set_yscale("log")
            img.axs[img_row][img_col].set_xscale("log")

            img.axs[img_row][img_col].set_xlim(0.8, 60)
            img.set_ticklabel_str(img_row, img_col, 1,[1,5,10,20,40], ["$1$","$5$","$10$","$20$","$40$"])

            # img.axs[img_row][img_col].set_xlim(0.8, 200)
            # img.set_ticklabel_str(img_row, img_col, 1,[1,5,10,100], ["$1$","$5$","$10$","$100$"])
            img.axs[img_row][img_col].set_ylim(5 * 10 ** (-7), 5 * 10 ** (-4))

            # if tag not in [0, 6, 11, 15, 18, 20]:
            #     # img.del_axis(img_row,img_col,[1])
            #     img.axs[img_row][img_col].set_xticklabels([])
            if img_col != 0:
                img.axs[img_row][img_col].set_yticklabels([])
                # img.axs[img_row][img_col].set_yticks([])

            tag += 1



img.save_img(pic_nm_m)
img.close_img()