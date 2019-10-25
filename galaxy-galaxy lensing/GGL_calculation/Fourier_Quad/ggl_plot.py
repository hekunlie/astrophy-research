import h5py
import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot

data_path = "D:/result/"

img = Image_Plot()
img.subplots(1,1)

ylabel = "$\Delta\Sigma \; [\\rm{h} \cdot \\rm{M_{\odot}} \cdot \\rm{pc^{-2}}]$"
xlabel = "$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$"

for i in [1,3]:
    h5f = h5py.File(data_path + "cmass_result_w_%d.hdf5"%i,"r")
    res = h5f["/result"].value
    dist = h5f["/mean_dist"].value[0]
    h5f.close()
    print(dist.shape,res.shape)
    img.axs[0][0].errorbar(dist,res[0], res[1], marker="s", c="C%d"%i, mfc="none",label="W-%d"%i, capsize=4)

img.axs[0][0].legend(fontsize=img.legend_size)
img.axs[0][0].set_yscale("log")
img.axs[0][0].set_xscale("log")
img.axs[0][0].set_ylim(0.01,200)
img.set_label(0,0,0,ylabel,size=img.xy_lb_size)
img.set_label(0,0,1,xlabel,size=img.xy_lb_size)
img.save_img(data_path+"result.png")
img.show_img()



# # plot the difference percentage
# img = Image_Plot()
# img.subplots(1,1)
#
# ylabel = "$\Delta\Sigma \; [\\rm{h} \cdot \\rm{M_{\odot}} \cdot \\rm{pc^{-2}}]$"
# xlabel = "$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$"
# data_path = "E:/works/Galaxy-Galaxy_lensing/CLUSTER/"
# dist_min, dist_max = 100, 0
# for i in range(1, 5):
#     data_path1 = "E:/works/Galaxy-Galaxy_lensing/CLUSTER/1/"
#     data_path2 = "E:/works/Galaxy-Galaxy_lensing/CLUSTER/3/"
#     h5f1 = h5py.File(data_path1 + "CFHT_cluster_result_w_%d.hdf5"%i,"r")
#     h5f2 = h5py.File(data_path2 + "CFHT_cluster_result_w_%d.hdf5"%i,"r")
#     res = (h5f1["/result"].value[0] - h5f2["/result"].value[0])/h5f2["/result"].value[0]
#     dist = h5f1["/mean_dist"].value[0]
#     if dist.min()<dist_min:
#         dist_min = dist.min()
#     if dist.max() > dist_max:
#         dist_max = dist.max()
#     h5f1.close()
#     h5f2.close()
#
#     img.axs[0][0].plot(dist,res, marker="s", c="C%d"%i, mfc="none",label="W-%d"%i)
# dx = (dist_max - dist_min)*0.05
# xs = [dist_min - dx, dist_max + dx*2]
# img.axs[0][0].plot([xs[0], xs[1]], [0,0], c="grey", alpha=0.8, linestyle="--")
# img.axs[0][0].plot([xs[0], xs[1]], [0.1,0.1], c="grey", alpha=0.8, linestyle="--")
# img.axs[0][0].plot([xs[0], xs[1]], [-0.1,-0.1], c="grey", alpha=0.8, linestyle="--")
# img.axs[0][0].plot([xs[0], xs[1]], [0.2,0.2], c="grey", alpha=0.8, linestyle="--")
# img.axs[0][0].plot([xs[0], xs[1]], [-0.2,-0.2], c="grey", alpha=0.8, linestyle="--")
# img.axs[0][0].legend(fontsize=img.legend_size)
#
# img.axs[0][0].set_xscale("log")
#
# # img.axs[0][0].set_title("$SNR \geq 3.5$",fontsize=img.xy_lb_size)
#
# img.set_label(0,0,0,"$(\Delta\Sigma_{snr \geq 3.5}-\Delta\Sigma_{snr \geq 5.5})/\Delta\Sigma_{snr \geq 5.5}$",size=img.xy_lb_size)
# img.set_label(0,0,1,xlabel,size=img.xy_lb_size)
# img.save_img(data_path+"cutoff_1v3.png")
# img.show_img()