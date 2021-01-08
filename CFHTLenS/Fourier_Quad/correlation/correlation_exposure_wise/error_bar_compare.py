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
theta_bin_num = 7#int(argv[2])
resample_num = 300#int(argv[3])
pts_num = int(theta_bin_num*(zbin_num**2+zbin_num)/2)
data_path = "E:/works/correlation/CFHT"

pic_nm = data_path + "/error_bar_compare.png"

img = Image_Plot(fig_x=6,fig_y=4)
img.subplots(1, 2)
img.set_style()

result_npz = data_path + "/result_cache_100.npz"
npz = numpy.load(result_npz)
xi_p_sig_bench,xi_m_sig_bench = npz["arr_3"][0],npz["arr_4"][0]
print(xi_p_sig_bench.shape)
print(xi_m_sig_bench.shape)
for tag, resample_num in enumerate([200,300]):
    result_npz = data_path + "/result_cache_%d.npz"%resample_num
    npz = numpy.load(result_npz)
    xi_p_sig,xi_m_sig = npz["arr_3"][0]/xi_p_sig_bench,npz["arr_4"][0]/xi_m_sig_bench
    x = numpy.arange(0,theta_bin_num)

    for i in range(int((zbin_num**2+zbin_num)/2)):
        lb1 = None
        lb2 = None
        if i == 0:
            lb1 = "$\sigma_{\\xi_{+}}$ Jack %d"%resample_num
            lb2 = "$\sigma_{\\xi_{-}}$ Jack %d"%resample_num
        img.axs[0][0].scatter(x+i*30, xi_p_sig[i*theta_bin_num:(i+1)*theta_bin_num],c="C%d"%tag,label=lb1)
        img.axs[0][1].scatter(x+i*30, xi_m_sig[i*theta_bin_num:(i+1)*theta_bin_num],c="C%d"%tag,label=lb2)
    # img.axs[0][1].set_title()

img.axs[0][0].legend()
img.axs[0][1].legend()
img.save_img(pic_nm)
img.show_img()