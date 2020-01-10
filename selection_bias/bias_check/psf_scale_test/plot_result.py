import matplotlib
matplotlib.use("Agg")
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import numpy
import h5py
import os

parent_path = "/mnt/ddnfs/data_users/hkli/bias_check/data_from_pi/psf_scale_test"

img = Image_Plot(fig_x=7,fig_y=2,ypad=0,xpad=0.2)
img.subplots(4,1)

m_scale = 100
c_scale = 10000
x_tick_pos = []
for i in range(7):
    noise_free_path = parent_path + "/data_%d/shear_result_noise_free_epsf_ori_pdf.hdf5"%i
    noisy_path = parent_path + "/data_%d/shear_result_noisy_cpp_epsf_ori_pdf.hdf5"%i
    print(i,os.path.exists(noise_free_path))
    print(os.path.exists(noisy_path))
    #
    # continue
    h5f = h5py.File(noise_free_path,"r")
    sym_mc_nf = h5f["/sym_mc"][()]
    mean_mc_nf = h5f["/mean_mc"][()]
    h5f.close()

    h5f = h5py.File(noisy_path,"r")
    sym_mc_n = h5f["/sym_mc"][()]
    mean_mc_n = h5f["/mean_mc"][()]
    h5f.close()
    print(sym_mc_nf)
    print(mean_mc_nf)
    print(sym_mc_n)
    print(mean_mc_n)
    x = i*2
    x_tick_pos.append(x)
    # noise-free
    # m1/2
    img.axs[0][0].errorbar(x, sym_mc_nf[0,0]*m_scale, sym_mc_nf[0,1]*m_scale,c="C1",marker="o")
    img.axs[0][0].errorbar(x+0.2, sym_mc_nf[1,0]*m_scale, sym_mc_nf[1,1]*m_scale,c="C1",marker="v")

    img.axs[0][0].errorbar(x+0.5, mean_mc_nf[0,0]*m_scale, mean_mc_nf[0,1]*m_scale,c="C2",marker="o")
    img.axs[0][0].errorbar(x+0.7, mean_mc_nf[1,0]*m_scale, mean_mc_nf[1,1]*m_scale,c="C2",marker="v")

    # c1/2
    img.axs[1][0].errorbar(x, sym_mc_nf[0,2]*c_scale, sym_mc_nf[0,3]*c_scale,c="C1",marker="o")
    img.axs[1][0].errorbar(x+0.2, sym_mc_nf[1,2]*c_scale, sym_mc_nf[1,3]*c_scale,c="C1",marker="v")

    img.axs[1][0].errorbar(x+0.5, mean_mc_nf[0,2]*c_scale, mean_mc_nf[0,3]*c_scale,c="C2",marker="o")
    img.axs[1][0].errorbar(x+0.7, mean_mc_nf[1,2]*c_scale, mean_mc_nf[1,3]*c_scale,c="C2",marker="v")

    # noisy
    # m1/2
    img.axs[2][0].errorbar(x, sym_mc_n[0,0]*m_scale, sym_mc_n[0,1]*m_scale,c="C1",marker="o")
    img.axs[2][0].errorbar(x+0.2, sym_mc_n[1,0]*m_scale, sym_mc_n[1,1]*m_scale,c="C1",marker="v")

    img.axs[2][0].errorbar(x+0.5, mean_mc_n[0,0]*m_scale, mean_mc_n[0,1]*m_scale,c="C2",marker="o")
    img.axs[2][0].errorbar(x+0.7, mean_mc_n[1,0]*m_scale, mean_mc_n[1,1]*m_scale,c="C2",marker="v")

    # c1/2
    img.axs[3][0].errorbar(x, sym_mc_n[0,2]*c_scale, sym_mc_n[0,3]*c_scale,c="C1",marker="o")
    img.axs[3][0].errorbar(x+0.2, sym_mc_n[1,2]*c_scale, sym_mc_n[1,3]*c_scale,c="C1",marker="v")

    img.axs[3][0].errorbar(x+0.5, mean_mc_n[0,2]*c_scale, mean_mc_n[0,3]*c_scale,c="C2",marker="o")
    img.axs[3][0].errorbar(x+0.7, mean_mc_n[1,2]*c_scale, mean_mc_n[1,3]*c_scale,c="C2",marker="v")

for i in range(2):
    xs = img.axs[2*i][0].set_xlim()
    img.axs[2*i][0].plot([xs[0],xs[1]],[0.3,0.3],ls="--",alpha=0.4,c="k")
    img.axs[2*i][0].plot([xs[0],xs[1]],[0,0],ls="--",alpha=0.4,c="k")
    img.axs[2*i][0].plot([xs[0],xs[1]],[-0.3,-0.3],ls="--",alpha=0.4,c="k")

    xs = img.axs[2*i+1][0].set_xlim()
    img.axs[2*i+1][0].plot([xs[0],xs[1]],[1,1],ls="--",alpha=0.4,c="k")
    img.axs[2*i+1][0].plot([xs[0],xs[1]],[0,0],ls="--",alpha=0.4,c="k")
    img.axs[2*i+1][0].plot([xs[0],xs[1]],[-1,-1],ls="--",alpha=0.4,c="k")

labels = ["noise_free\n$10^2m_{1/2}$","noise_free\n$10^4c_{1/2}$",
          "noisy\n$10^2m_{1/2}$","noisy\n$10^4c_{1/2}$"]
for i in range(4):
    img.set_label(i,0,0,labels[i],fontsize=img.xy_lb_size-2)

psf_scale = ["$r_d=1.6$", "$r_d=2$", "$r_d=2.5$",
             "$r_d=3$", "$r_d=4$", "$r_d=6$", "$r_d=7$"]
# img.axs[3][0].set_xticks(x_tick_pos)
# img.axs[3][0].set_xticklabels(psf_scale)
img.set_ticklabel_str(3,0,1,x_tick_pos,psf_scale)
img.save_img(parent_path+"/result.png")

