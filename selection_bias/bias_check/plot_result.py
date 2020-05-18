import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot


c_scale = 1000
m_scale = 100

names = ["%d"%i for i in range(10)]
# names = []
# names.append("total")
data_type = "noise_free"

for nm in names:
    parent_path = "E:/data/weight_test/step_1_psf_4_large_sample_48x48"
    data_path = parent_path + "/%s"%nm
    pic_nm = parent_path + "/%s_%s.png"%(data_type,nm)
    img = Image_Plot(xpad=0.2,ypad=0.1,fig_x=5,fig_y=4)
    img.subplots(2,3)

    # weight_label = [5,6,7]
    # weight_pow = [1,2,3]
    weight_label = [4,5,6]
    weight_pow = [1,2,3,4]

    xlabels = ["$w^%d$"%i for i in range(len(weight_pow)+1)]
    # weight = ["$w=\\frac{1}{P_{k0,fit}}$",
    #           "$w=\\frac{1}{Max(P_{k0},P_{k0,fit})}$", "$w=\\frac{1}{Flux_{true}}$"]
    weight = ["$w=\\frac{1}{P_{k0}}$","$w=\\frac{1}{P_{k0,fit}}$",
              "$w=\\frac{1}{Max(P_{k0},P_{k0,fit})}$"]
    print(xlabels)

    for i,w in enumerate(weight_label):
        ave_mcs = numpy.zeros((8, len(weight_pow)+1))
        pdf_mcs = numpy.zeros((8, len(weight_pow)+1))

        x = numpy.arange(0,len(weight_pow)+1)
        h5f = h5py.File(data_path + "/shear_result_%s_epsf_0_flux_0.hdf5"%data_type,"r")
        mean_mc = h5f["/mean_mc"][()].flatten()
        pdf_mc = h5f["/sym_mc"][()].flatten()
        h5f.close()

        ave_mcs[:,0] = mean_mc
        pdf_mcs[:,0] = pdf_mc

        for j, n in enumerate(weight_pow):
            h5f = h5py.File(data_path + "/shear_result_%s_epsf_%d_flux_%d.hdf5"%(data_type,w,n), "r")
            mean_mc = h5f["/mean_mc"][()].flatten()
            pdf_mc = h5f["/sym_mc"][()].flatten()
            h5f.close()

            ave_mcs[:,1+j] = mean_mc
            pdf_mcs[:,1+j] = pdf_mc


        img.axs[0][i].errorbar(x, ave_mcs[0]*m_scale,ave_mcs[1]*m_scale, label="Ave. $10^2m1$",
                               capsize=img.cap_size, marker="s")
        img.axs[0][i].errorbar(x, ave_mcs[4]*m_scale,ave_mcs[5]*m_scale, label="Ave. $10^2m2$",
                               capsize=img.cap_size, marker="s")
        img.axs[0][i].errorbar(x+0.15, pdf_mcs[0]*m_scale,pdf_mcs[1]*m_scale, label="PDF $10^2m1$",
                               capsize=img.cap_size, marker="o",ls="--")
        img.axs[0][i].errorbar(x+0.15, pdf_mcs[4]*m_scale,pdf_mcs[5]*m_scale, label="PDF $10^2m2$",
                               capsize=img.cap_size, marker="o",ls="--")

        img.axs[1][i].errorbar(x, ave_mcs[2]*c_scale,ave_mcs[3]*c_scale, label="Ave. $10^4c1$",
                               capsize=img.cap_size, marker="s")
        img.axs[1][i].errorbar(x, ave_mcs[6]*c_scale,ave_mcs[7]*c_scale, label="Ave. $10^4c2$",
                               capsize=img.cap_size, marker="s")
        img.axs[1][i].errorbar(x+0.15, pdf_mcs[2]*c_scale,pdf_mcs[3]*c_scale, label="PDF $10^4c1$",
                               capsize=img.cap_size, marker="o",ls="--")
        img.axs[1][i].errorbar(x+0.15, pdf_mcs[6]*c_scale,pdf_mcs[7]*c_scale, label="PDF $10^4c2$",
                               capsize=img.cap_size, marker="o",ls="--")

        img.axs[0][i].legend(ncol=2,fontsize=img.legend_size-8,loc="upper left")
        img.axs[1][i].legend(ncol=2,fontsize=img.legend_size-8,loc="upper left")
        img.set_ticklabel_str(0,i,1,x,xlabels)
        img.set_ticklabel_str(1,i,1,x,xlabels)
        img.axs_text(0,i,0.1,0.1,weight[i],text_color="k")
        img.axs_text(1,i,0.1,0.1,weight[i],text_color="k")
        img.axs[0][i].set_ylim(-1.5,1.5)
        img.axs[1][i].set_ylim(-0.5,0.5)
    img.save_img(pic_nm)
    # img.show_img()
    img.close_img()