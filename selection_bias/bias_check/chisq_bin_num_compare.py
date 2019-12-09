import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box



def mc_plot(result_arr, mc_arr, img_name):

    img = Image_Plot(fig_x=6, fig_y=4, xpad=0.25, ypad=0.2)
    img.subplots(1, 2)

    scale = 1000

    img.axs[0][0].errorbar(result_arr[0], result_arr[1], result_arr[2],
                           marker="s", fmt=" ", label="g1", capsize=img.cap_size)
    img.axs[0][0].errorbar(result_arr[3], result_arr[4], result_arr[5],
                           marker="s", fmt=" ", label="g2", capsize=img.cap_size)
    img.axs[0][1].errorbar(result_arr[0], scale*(result_arr[1]-result_arr[0]), scale*result_arr[2],
                           marker="s", fmt=" ", label="g1", capsize=img.cap_size)
    img.axs[0][1].errorbar(result_arr[3], scale*(result_arr[4]-result_arr[3]), scale*result_arr[5],
                           marker="s", fmt=" ", label="g2", capsize=img.cap_size)

    mc_str = "Mean\n$m_1 = %.5f (%.5f)$\n$c_1 = %.5f (%.5f)$\n$m_2 = %.5f (%.5f)$\n$c_2 = %.5f (%.5f)$" \
             % (mc_arr[0, 0], mc_arr[0, 1], mc_arr[0, 2], mc_arr[0, 3],
                mc_arr[1, 0], mc_arr[1, 1], mc_arr[1, 2], mc_arr[1, 3])
    img.axs_text(0, 0, 0.8, 0.05, mc_str, text_fontsize=img.legend_size)

    img.set_label(0, 0, 0, "EST g")
    img.set_label(0, 1, 1, "$10^3 (\hat{g} - g_{true})$")
    img.axs[0][0].plot([-0.1, 0.1], [-0.1, 0.1], ls="dashed", c="grey", alpha=0.5)
    img.axs[0][0].set_ylim(-0.05, 0.05)
    img.axs[0][1].set_ylim(-0.1, 0.1)
    for i in range(2):
        img.set_label(0, i, 1, "TRUE g")
        img.axs[0][i].legend(fontsize=img.legend_size, loc="lower right")
        img.axs[0][i].set_xlim(-0.05, 0.05)
    img.save_img(img_name)
    img.close_img()

# bin_num = [2,8,12,16,20,40,80]
fit_range =[0.005, 0.01, 0.015, 0.02, 0.025, 0.03, 0.035, 0.04, 0.045]
bin_num = [2,4,8,12,16,20,40,60,80,120,160,200]

fit_id = 0
for fit_id in range(9):
    fit_label = "fit_%d"%fit_id
    parent_path = "E:/works/Group_meeting/2019-11-25-shear_bias_checking/result/data_from_pi/3/result"
    mc1 = numpy.zeros((len(bin_num), 4))
    mc2 = numpy.zeros((len(bin_num), 4))
    for ig, num in enumerate(bin_num):
        h5f = h5py.File(parent_path + "/result_noisy_bin_num_%d_%s.hdf5"%(num,fit_label),"r")
        sym_result_array = h5f["/sym_result"].value
        sym_mc_array = h5f["/sym_mc"].value
        h5f.close()
        mc1[ig] = sym_mc_array[0]
        mc2[ig] = sym_mc_array[1]

        img1_path = parent_path + "/result_noisy_bin_num_%d_%s.png"%(num,fit_label)
        mc_plot(sym_result_array,sym_mc_array,img1_path)


    img = Image_Plot(fig_x=8, fig_y=6,xpad=0.25,ypad=0.2)
    img.subplots(1,2)
    scale = [100,10000]
    lb = [["m1","m2"],["c1","c2"]]
    ylabels = ["$10^2m$","$10^4c$"]
    for i in range(2):
        img.axs[0][i].errorbar(range(len(bin_num)),mc1[:,i*2]*scale[i], mc1[:,2*i+1]*scale[i],
                               label=lb[i][0],capsize=img.cap_size,marker="s")
        img.axs[0][i].errorbar(range(len(bin_num)),mc2[:,i*2]*scale[i], mc2[:,2*i+1]*scale[i],
                               label=lb[i][1],capsize=img.cap_size,marker="s")
        img.set_label(0,i,0,ylabels[i])
        img.axs[0][i].legend(fontsize=img.legend_size)
        img.axs[0][i].set_title("fit_range: g$\pm$%.4f"%fit_range[fit_id])
        xs = img.axs[0][i].set_xlim()
        img.axs[0][i].plot([xs[0],xs[1]],[0,0],ls="--",alpha=0.5,c="grey")
        img.axs[0][i].plot([xs[0],xs[1]],[-0.1,-0.1],ls="--",alpha=0.5,c="grey")
        img.axs[0][i].plot([xs[0],xs[1]],[-0.2,-0.2],ls="--",alpha=0.5,c="grey")
        img.axs[0][i].plot([xs[0],xs[1]],[-0.3,-0.3],ls="--",alpha=0.5,c="grey")

    y1 = mc1[:,0] - mc1[:,1]*1.5
    y2 = mc2[:,0] - mc2[:,1]*1.5
    for i in range(len(bin_num)):
        img.axs[0][0].text(i,min(y1[i],y2[i])*scale[0],"%d bins"%bin_num[i],color="green", ha="left",  va="center",
                           rotation=60,fontsize=img.legend_size-4)
    for i in range(2):
        img.del_tick(0,i,[1])
    img.save_img(parent_path+"/result_noisy_bin_num_compare_%s.png"%fit_label)
    img.close_img()
    # img.show_img()