import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
# path.append('%s/work/mylib/' % my_home)
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import h5py
import numpy

zbin_num = 6#int(argv[1])
theta_bin_num = 5#int(argv[2])
resample_num = 300#int(argv[3])


h5f = h5py.File("D:/result.hdf5","r")
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
    xi_p_sig[:,i] = xi_p_sub[:,i].std()
    xi_m_sig[:,i] = xi_m_sub[:,i].std()
    img = Image_Plot(fig_x=6, fig_y=4)
    img.subplots(1, 1)
    # img.axs[0][0].hist(xi_p_sub[:,i],50)
    img.axs[0][0].scatter(range(resample_num),xi_p_sub[:,i])
    img.axs[0][0].set_yscale("symlog")
    # img.axs[0][0].set_ylim(10**(-7),10**(-4))
    img.show_img()

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
img = Image_Plot(fig_x=3, fig_y=2,xpad=0,ypad=0)
img.subplots(zbin_num, zbin_num)


tag = 0
for i in range(zbin_num):
    for j in range(zbin_num):
        img_row = zbin_num - j - 1
        img_col = i
        if j >= i:
            st, ed = int(tag*theta_bin_num), int((tag+1)*theta_bin_num)
            # print(theta[st:ed])

            img.axs[img_row][img_col].errorbar(theta[0,st:ed], xi_p[0,st:ed],yerr=xi_p_sig[0,st:ed],marker="s",ms=5,label="$\\xi_{+}$",capsize=3)
            img.axs[img_row][img_col].errorbar(theta[0,st:ed], xi_m[0,st:ed],yerr=xi_m_sig[0,st:ed],marker="s",ms=5,label="$\\xi_{-}$",capsize=3)
            img.axs[img_row][img_col].legend(loc="upper right",fontsize=18)
            img.axs[img_row][img_col].set_yscale("log")
            img.axs[img_row][img_col].set_xscale("log")

            if tag not in [0,6,11,15,18,20]:
                img.del_axis(img_row,img_col,[1])
            if img_col != 0:
                img.del_ticklabel(img_row,img_col,[0])

            img.axs_text(img_row, img_col, 0.8, 0.1, "%d-%d"%(i+1,j+1),text_fontsize=18,text_color="k")
            img.axs[img_row][img_col].set_ylim(10**(-7),5*10**(-4))
            tag += 1
    for j in range(zbin_num-i,zbin_num):
        img.figure.delaxes(img.axs[i][j])
img.save_img("D:/result.png")
img.close_img()
# img.show_img()
