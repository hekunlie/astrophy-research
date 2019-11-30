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


mode = "chisq"

fq = Fourier_Quad(12,124)
g1 = numpy.array([-0.04, 0, 0.04])
g2 = numpy.array([0.04, -0.04, 0])
shear_num = g1.shape[0]

sub_num = 4
data_col = 4
result_g1 = numpy.zeros((5, shear_num))
result_g1[4] = g1
result_g2 = numpy.zeros((5, shear_num))
result_g2[4] = g2

parent_path = "E:/test/data_1"

if mode == "all":
    for ir in range(1):

        img = Image_Plot(fig_x=5, fig_y=4, cap_size=2.5)
        img.subplots(2, shear_num + 2)

        pic_nm = parent_path + "/data_total_noise-free.png"#%ir
        result_path = parent_path + "/data_total_noise-free.txt"#%ir
        for tag in range(shear_num):

            # h5f = h5py.File(parent_path + "/data_rotation_%d/data_%d.hdf5"%(ir,tag),"r")
            h5f = h5py.File(parent_path + "/data_total/data_%d.hdf5"%tag,"r")
            data = h5f["/data"].value
            print(data.shape)
            h5f.close()

            data_nms = ["Noise-free", "Noise_power", "Noise_residual", "Normal"]

            # img = Image_Plot(fig_x=5, fig_y=3)
            # img.subplots(1, 4)

            for i in range(1):
                left, right = -0.1, 0.1
                if i == 1:
                    left, right = -0.01, 0.01

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
                g1_sym, g1_sig_sym = fq.find_shear(mg1, mnu1, 8,left=left,right=right,fit_num=20, chi_gap=20,fig_ax=img.axs[0][tag])[:2]
                g2_sym, g2_sig_sym = fq.find_shear(mg2, mnu2, 8,left=left,right=right, fit_num=20, chi_gap=20,fig_ax=img.axs[1][tag])[:2]

                result_g1[0,tag] = g1_mean
                result_g1[1,tag] = g1_sig_mean
                result_g1[2,tag] = g1_sym
                result_g1[3,tag] = g1_sig_sym

                result_g2[0, tag] = g2_mean
                result_g2[1, tag] = g2_sig_mean
                result_g2[2, tag] = g2_sym
                result_g2[3, tag] = g2_sig_sym
                str1 = "%s\nTrue： g1: %.5f, g2: %.5f"%(data_nms[i],g1[tag], g2[tag])
                str2 = "\nMEAN: g1: %.5f (%.5f) g2: %.5f (%.5f)"%(g1_mean, g1_sig_mean, g2_mean, g2_sig_mean)
                str3 = "\nDiff: g1: %.5f (%.5f) g2: %.5f (%.5f)"%(g1[tag]-g1_mean, g1_sig_mean, g2[tag]-g2_mean, g2_sig_mean)
                str4 = "\nSYM: g1: %.5f (%.5f) g2: %.5f (%.5f)"%(g1_sym, g1_sig_sym,g2_sym, g2_sig_sym)
                str5 = "\nDiff: g1: %.5f (%.5f) g2: %.5f (%.5f)"%(g1[tag]-g1_sym, g1_sig_sym,g2[tag]-g2_sym, g2_sig_sym)
                text_str = str1+str2+str3+str4+str5
                print(text_str)

                # img.axs[1][0].hist(mg1, 100)
                # img.axs[1][1].hist(mg1, 100)
                # img.axs[0][2].hist(mn, 100)
                # img.axs[0][3].hist(mu, 100)
                # img.axs[1][2].hist(mnu1, 100)
                # img.axs[1][3].hist(mnu2, 100)

                #
                # result_g1[tag,i*sub_num:(i+1)*sub_num] = g1_mean, g1_sig_mean,g1_sym, g1_sig_sym
                # result_g2[tag,i*sub_num:(i+1)*sub_num] = g2_mean, g2_sig_mean,g2_sym, g2_sig_sym

            #     img.axs[0][i].hist(mg1, 1000)
            #     img.axs_text(0,i,0.9, 0.1,text_str)
            # img.show_img()
            # img.close_img()

        img.axs[0][shear_num].errorbar(g1, result_g1[0], result_g1[1], c="C1",fmt=" ",label="Mean g1",capsize=img.cap_size)
        img.axs[0][shear_num].scatter(g1, result_g1[0], c="C1")
        img.axs[0][shear_num].errorbar(g1, result_g1[2], result_g1[3], c="C2",fmt=" ",label="SYM g1",capsize=img.cap_size)
        img.axs[0][shear_num].scatter(g1, result_g1[2], c="C2")

        img.axs[1][shear_num].errorbar(g2, result_g2[0], result_g2[1],fmt=" ",label="Mean g2",capsize=img.cap_size)
        img.axs[1][shear_num].scatter(g2, result_g2[0],c="C1")
        img.axs[1][shear_num].errorbar(g2, result_g2[2], result_g2[3],fmt=" ",label="SYM g2",capsize=img.cap_size)
        img.axs[1][shear_num].scatter(g2, result_g2[2], c="C2")
        img.axs[0][shear_num].legend(fontsize=img.legend_size)
        img.axs[1][shear_num].legend(fontsize=img.legend_size)

        scale = 1000

        dg1_mean = numpy.abs(result_g1[0]) - numpy.abs(g1)
        dg1_sym = numpy.abs(result_g1[2]) - numpy.abs(g1)

        dg2_mean = numpy.abs(result_g2[0]) - numpy.abs(g2)
        dg2_sym = numpy.abs(result_g2[2]) - numpy.abs(g2)

        img.axs[0][shear_num+1].errorbar(g1, scale*dg1_mean, scale*result_g1[1],c="C1",
                                         fmt=" ",label="$10^3 |dg1|,Mean$",capsize=img.cap_size)
        img.axs[0][shear_num+1].scatter(g1, scale*dg1_mean,c="C1",)

        img.axs[0][shear_num+1].errorbar(g1, scale*dg1_sym, scale*result_g1[3],c="C2",
                                         fmt=" ",label="$10^3 |dg1|,SYM$",capsize=img.cap_size)
        img.axs[0][shear_num+1].scatter(g1, scale*dg1_sym,c="C2")

        img.axs[1][shear_num+1].errorbar(g2, scale*dg2_mean, scale*result_g2[1],c="C1",
                                         fmt=" ",label="$10^3 |dg2|,Mean$",capsize=img.cap_size)
        img.axs[1][shear_num+1].scatter(g2, scale*dg2_mean,c="C1")

        img.axs[1][shear_num+1].errorbar(g2, scale*dg2_sym, scale*result_g2[3],c="C2",
                                         fmt=" ",label="$10^3 |dg2|,SYM$",capsize=img.cap_size)
        img.axs[1][shear_num+1].scatter(g2, scale*dg2_sym,c="C2")

        img.axs[0][shear_num+1].legend(fontsize=img.legend_size)
        img.axs[1][shear_num+1].legend(fontsize=img.legend_size)

        img.save_img(pic_nm)
        # img.show_img()
        img.close_img()
        header = "row: g1-mean, dg1-mean, g1-pdf, dg1-pdf, g1-true, g2-mean, dg2-mean, g2-pdf, dg2-pdf, g2-true"
        result = numpy.row_stack((result_g1, result_g2))
        numpy.savetxt(result_path, X=result, header=header)
        print(ir)

else:

    shear_tag = 0

    pic_nm = parent_path + "/data_total_chisq_%d.png"%shear_tag
    result_path = parent_path + "/data_total_chisq_%d.txt"%shear_tag

    h5f = h5py.File(parent_path + "/data_total/data_%d.hdf5" % shear_tag, "r")
    data = h5f["/data"].value
    print(data.shape)
    h5f.close()

    labels = ["Noise-free", "Noise_image", "Noise_power_residual", "Noisy_image"]
    img = Image_Plot(fig_x=5, fig_y=3, cap_size=2.5,xpad=0.2,ypad=0.25)
    img.subplots(2, 2)


    for data_tag in [0,2,3]:

        mg1 = data[:, data_tag*sub_num]
        mg2 = data[:, data_tag*sub_num+1]
        mn = data[:, data_tag*sub_num+2]
        mu = data[:, data_tag*sub_num+3]
        mnu1 = mn + mu
        mnu2 = mn - mu

        dg = 0.0015
        if data_tag == 1:
            g1_guess = numpy.linspace(-dg,  dg, 21)
            g2_guess = numpy.linspace(-dg,  dg, 21)
        elif data_tag == 2:
            g1_guess = numpy.linspace(-0.1, 0.1, 21)
            g2_guess = numpy.linspace(-0.1, 0.1, 21)
        else:
            g1_guess = numpy.linspace(g1[shear_tag] - dg, g1[shear_tag] + dg, 21)
            g2_guess = numpy.linspace(g2[shear_tag] - dg, g2[shear_tag] + dg, 21)

        g1_hat, g1_chisq = fq.get_chisq_range(mg1, mnu1, 8, g1_guess, limited=40)
        g2_hat, g2_chisq = fq.get_chisq_range(mg2, mnu2, 8, g2_guess, limited=40)

        img.axs[0][0].plot(g1_hat, g1_chisq, marker="o", c="C%d"%data_tag, label=labels[data_tag])
        img.axs[1][0].plot(g2_hat, g2_chisq, marker="o", c="C%d"%data_tag, label=labels[data_tag])

        if data_tag != 2:
            # estimate shear signal
            coeff = tool_box.fit_1d(g1_hat, g1_chisq, 2, "scipy")
            g1_estimate = -coeff[1] / 2. / coeff[2]
            g1_estimate_sig = 0.70710678118 / numpy.sqrt(coeff[2])
            # estimate shear signal
            coeff = tool_box.fit_1d(g2_hat, g2_chisq, 2, "scipy")
            g2_estimate = -coeff[1] / 2. / coeff[2]
            g2_estimate_sig = 0.70710678118 / numpy.sqrt(coeff[2])

            # plot chi squared
            img.axs[0][1].plot(g1_hat, g1_chisq, marker="o", c="C%d"%data_tag, label=labels[data_tag])
            ys = img.axs[0][1].set_ylim()
            # 1 sigma range
            img.axs[0][1].fill_between([g1_estimate-g1_estimate_sig, g1_estimate+g1_estimate_sig],
                                       ys[0],ys[1], facecolor="C%d"%data_tag, alpha=0.2)
            # true signal
            img.axs[0][1].plot([g1[shear_tag],g1[shear_tag]],[ys[0],ys[1]],ls="dashed",c="grey")
            img.axs[0][1].plot([g1_estimate, g1_estimate],[ys[0],ys[1]],ls="dotted",c="C%d"%data_tag)

            img.axs[1][1].plot(g2_hat, g2_chisq, marker="o", c="C%d"%data_tag, label=labels[data_tag])
            ys = img.axs[1][1].set_ylim()
            img.axs[1][1].fill_between([g2_estimate-g2_estimate_sig, g2_estimate+g2_estimate_sig],
                                       ys[0],ys[1], facecolor="C%d"%data_tag,alpha=0.2)
            img.axs[1][1].plot([g2[shear_tag],g2[shear_tag]],[ys[0],ys[1]],ls="dashed",c="grey")
            img.axs[1][1].plot([g2_estimate, g2_estimate],[ys[0],ys[1]],ls="dotted",c="C%d"%data_tag)

    for i in range(4):
        m,n = divmod(i,2)
        img.set_label(m,n,0,"$\chi^2$")
        img.set_label(m,n,1,"$\hat{g}$")
        img.axs[m][n].legend(fontsize=img.legend_size)
    img.save_img(pic_nm)
    img.show_img()
