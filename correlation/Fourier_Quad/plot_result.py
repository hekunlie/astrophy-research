import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
# path.append('%s/work/mylib/' % my_home)
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import h5py
import numpy
import matplotlib.pyplot as plt


zbin_num = 6#int(argv[1])
theta_bin_num = 5#int(argv[2])
resample_num = 300#int(argv[3])

pic_nm = "D:/result_%d.png"%resample_num
h5f = h5py.File("D:/result_%d.hdf5"%resample_num,"r")
print(list(h5f.keys()))
print(list(h5f["/0"].keys()))

xi_p = (-h5f["/0/tt"][()] - h5f["/0/xx"][()]).reshape((1,105))
xi_m = (-h5f["/0/tt"][()] + h5f["/0/xx"][()]).reshape((1,105))
theta = h5f["/0/theta"][()].reshape((1,105))

xi_p_sig = numpy.zeros_like(xi_p)
xi_m_sig = numpy.zeros_like(xi_p)

xi_p_sub = numpy.zeros((resample_num, 105))
xi_m_sub = numpy.zeros((resample_num, 105))

for i in range(1,resample_num+1):

    xi_p_sub[i-1] = (-h5f["/%d/tt"%i][()] - h5f["/%d/xx"%i][()]).reshape((1,105))
    xi_m_sub[i-1] = (-h5f["/%d/tt"%i][()] + h5f["/%d/xx"%i][()]).reshape((1,105))

for i in range(105):
    xi_p_sig[:,i] = xi_p_sub[:,i].std()*numpy.sqrt(resample_num-1)
    xi_m_sig[:,i] = xi_m_sub[:,i].std()*numpy.sqrt(resample_num-1)
    # img = Image_Plot(fig_x=6, fig_y=4)
    # img.subplots(1, 1)
    # # img.axs[0][0].hist(xi_p_sub[:,i],50)
    # img.axs[0][0].scatter(range(resample_num),xi_p_sub[:,i])
    # img.axs[0][0].set_yscale("symlog")
    # # img.axs[0][0].set_ylim(10**(-7),10**(-4))
    # img.show_img()
# print(xi_p_sig)
# y = -h5f["/0/xx"][()].reshape((1,105))[0,:5]
# print(y.max(), y.min())
# img = Image_Plot(fig_x=6, fig_y=4)
# img.subplots(1,1)
# img.axs[0][0].scatter(range(105)[:5],y)
# img.axs[0][0].set_yscale("log")
# img.axs[0][0].set_ylim(10**(-7),10**(-4))
# img.show_img()
# img.clos_img()
# h5f.close()
# exit()


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

            img.axs[img_row][img_col].errorbar(theta[0,st:ed], xi_p[0,st:ed],yerr=xi_p_sig[0,st:ed],
                                               lw=img.plt_line_width,elinewidth=img.plt_line_width,
                                               marker="o",fmt=" ",mfc="none",ms=12,label="$\\xi_{+}$",capsize=5)
            img.axs[img_row][img_col].errorbar(theta[0,st:ed], xi_m[0,st:ed],yerr=xi_m_sig[0,st:ed],
                                               lw=img.plt_line_width,elinewidth=img.plt_line_width,
                                               marker="o",fmt=" ",mfc="none",ms=12,label="$\\xi_{-}$",capsize=5)
            if tag == 0:
                img.axs[img_row][img_col].legend(loc="lower left",bbox_to_anchor=(1.2,0.25),
                                                 fancybox=False,fontsize=img.legend_size)

            img.axs_text(img_row, img_col, 0.8, 0.8, "%d-%d" % (i + 1, j + 1), text_fontsize=img.legend_size, text_color="k")

            # img.axs[img_row][img_col].set_xticks(fontsize=img.legend_size)
            # img.axs[img_row][img_col].set_yticks(fontsize=img.legend_size)
            img.axs[img_row][img_col].set_yscale("log")
            img.axs[img_row][img_col].set_xscale("log")

            img.axs[img_row][img_col].set_xlim(1.1, 60)
            img.axs[img_row][img_col].set_ylim(5 * 10 ** (-7), 5 * 10 ** (-4))

            # if tag not in [0, 6, 11, 15, 18, 20]:
            #     # img.del_axis(img_row,img_col,[1])
            #     img.axs[img_row][img_col].set_xticklabels([])
            if img_col != 0:
                img.axs[img_row][img_col].set_yticklabels([])
                # img.axs[img_row][img_col].set_yticks([])

            tag += 1



img.save_img(pic_nm)
img.close_img()
# img.show_img()
