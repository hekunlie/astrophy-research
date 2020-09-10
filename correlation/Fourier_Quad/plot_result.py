import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
# path.append('%s/work/mylib/' % my_home)
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import h5py
import numpy
import tool_box

zbin_num = 6#int(argv[1])
theta_bin_num = 5#int(argv[2])
resample_num = 200#int(argv[3])

img = Image_Plot(fig_x=4, fig_y=3,xpad=0.3,ypad=0.3)
img.subplots(zbin_num, zbin_num)

h5f = h5py.File("D:/result.hdf5","r")
chi_tt = h5f["/0/chi_tt"][()]
chi_tt_sig = h5f["/0/chi_tt_sig"][()]
chi_xx = h5f["/0/chi_xx"][()]
chi_xx_sig = h5f["/0/chi_xx_sig"][()]
theta = h5f["/0/theta"][()]

chi_p = chi_tt + chi_xx
chi_m = chi_tt - chi_xx
chi_sig = numpy.sqrt(chi_tt_sig**2 + chi_xx**2)
print(chi_tt.shape)
print(chi_xx.shape)
print(theta.shape)

theta = tool_box.set_bin_log(0.8, 60, theta_bin_num+1)
theta = (theta[1:] + theta[:-1])/2
tag = 0
for i in range(zbin_num):
    for j in range(zbin_num):
        img_row = zbin_num - j - 1
        img_col = i
        if j >= i:
            st, ed = int(tag*theta_bin_num), int((tag+1)*theta_bin_num)
            # print(theta[st:ed])

            img.axs[img_row][img_col].errorbar(theta, chi_p[st:ed],yerr=chi_sig[st:ed],label="$\chi_{+}$",capsize=3)
            img.axs[img_row][img_col].errorbar(theta, chi_m[st:ed],yerr=chi_sig[st:ed],label="$\chi_{-}$",capsize=3)
            img.axs[img_row][img_col].legend(loc="upper right")
            img.axs[img_row][img_col].set_yscale("log")
            img.axs[img_row][img_col].set_xscale("log")
            # if tag not in [0,6,11,15,18,20]:
            #     img.del_axis(img_row,img_col,[1])
            img.axs_text(img_row, img_col, 0.8, 0.1, "%d-%d"%(i+1,j+1),text_fontsize=18,text_color="k")
            tag += 1
        else:
            img.del_axis(img_row,img_col,[0,1],[0,1,2,3])
img.save_img("D:/result.png")
img.show_img()
h5f.close()