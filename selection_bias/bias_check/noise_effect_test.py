# import matplotlib
# matplotlib.use("Agg")
import numpy
from sys import path
# path.append('%s/work/fourier_quad/'%my_home)
path.append("D:/Github/astrophy-research/mylib/")
import tool_box
from plot_tool import Image_Plot
import h5py
from Fourier_Quad import Fourier_Quad


fq = Fourier_Quad(12,124)
g1 = numpy.array([-0.04, 0, 0.04])
g2 = numpy.array([0.04, -0.04, 0])
shear_num = g1.shape[0]

sub_num = 4
data_col = 4
result_g1 = numpy.zeros((shear_num, data_col*sub_num))
result_g2 = numpy.zeros((shear_num, data_col*sub_num))

for tag in range(2,3):#shear_num):
    h5f = h5py.File("D:\\result\data\\data_%d.hdf5"%tag,"r")
    data = h5f["/data"].value
    print(data.shape)
    h5f.close()

    # data_noise_free = data[:,:4]
    # data_noise_power = data[:,4:8]
    # data_noise_residual = data[:,8:12]
    # data_normal = data[:,12:]
    # datas = [data_noise_free, data_noise_power, data_noise_residual, data_normal]
    data_nms = ["Noise-free", "Noise_power", "Noise_residual", "Normal"]

    # img = Image_Plot(fig_x=5, fig_y=3)
    # img.subplots(1, 4)

    for i in range(sub_num):
        left, right = -0.1, 0.1
        if i == 1:
            left, right = -0.01,0.01

        img = Image_Plot(fig_x=5, fig_y=3)
        img.subplots(2, 4)

        mg1 = data[:, sub_num*i]
        mg2 = data[:, sub_num*i + 1]
        mn = data[:, sub_num*i + 2]
        mu = data[:, sub_num*i + 3]
        mnu1 = mn + mu
        mnu2 = mn - mu

        # img.show_img()
        # img.close_img()
        g1_mean, g1_sig_mean = fq.find_shear_mean(mg1, mn)
        g2_mean, g2_sig_mean = fq.find_shear_mean(mg2, mn)
        g1_sym, g1_sig_sym = fq.find_shear(mg1, mnu1, 6,left=left,right=right,fit_num=20, chi_gap=20,fig_ax=img.axs[0][0])[:2]
        g2_sym, g2_sig_sym = fq.find_shear(mg2, mnu2, 6,left=left,right=right, fit_num=20, chi_gap=20,fig_ax=img.axs[0][1])[:2]

        text_str = "%s\nTrue： g1: %.5f, g2: %.5f\nMEAN: g1: %.5f (%.5f) g2: %.5f (%.5f)\nDiff: g1: %.5f (%.5f) g2: %.5f (%.5f)" \
                   "\nSYM: g1: %.5f (%.5f) g2: %.5f (%.5f)\nDiff: g1: %.5f (%.5f) g2: %.5f (%.5f)"\
                   %(data_nms[i],g1[tag], g2[tag],g1_mean, g1_sig_mean, g2_mean, g2_sig_mean, g1_mean-g1[tag], g1_sig_mean, g2_mean-g2[tag], g2_sig_mean,
                     g1_sym, g1_sig_sym,g2_sym, g2_sig_sym, g1_sym-g1[tag], g1_sig_sym,g2_sym-g2[tag], g2_sig_sym)
        print(text_str)

        img.axs[1][0].hist(mg1, 100)
        img.axs[1][1].hist(mg1, 100)
        img.axs[0][2].hist(mn, 100)
        img.axs[0][3].hist(mu, 100)
        img.axs[1][2].hist(mnu1, 100)
        img.axs[1][3].hist(mnu2, 100)
        img.show_img()
        img.close_img()
        #
        # result_g1[tag,i*sub_num:(i+1)*sub_num] = g1_mean, g1_sig_mean,g1_sym, g1_sig_sym
        # result_g2[tag,i*sub_num:(i+1)*sub_num] = g2_mean, g2_sig_mean,g2_sym, g2_sig_sym

    #     img.axs[0][i].hist(mg1, 1000)
    #     img.axs_text(0,i,0.9, 0.1,text_str)
    # img.show_img()
    # img.close_img()
