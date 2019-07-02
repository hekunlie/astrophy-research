import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import matplotlib.pyplot as plt
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
import h5py


fore_source = argv[1]
data_file = argv[2]

total_path = "/mnt/perc/hklee/CFHT/gg_lensing/"

result_path = total_path + "result/%s/fourier/result/fourier_%s_result_%s.hdf5" % (fore_source, fore_source, data_file)
dens_pic_path = total_path + "result/%s/fourier/result/fourier_%s_%s" % (fore_source, fore_source, data_file)
dens_r_pic_path = total_path + "result/%s/fourier/result/fourier_%s_%s_sigmaxr" % (fore_source, fore_source, data_file)

h5f = h5py.File(result_path, "r")
result = h5f["/data"].value
h5f.close()


ylabels = ["$\gamma$", "$\Delta\Sigma \; [\\rm{M_{\odot}} \cdot \\rm{pc^{-2}}]$"]
ylabels_r = "$\\rm{R}\Delta\Sigma \; [\\rm{10^6\ M_{\odot}} \cdot \\rm{pc^{-2}}]$"
xlabel = "$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$"

gt_lb, gx_lb = 0, 2
sigt_lb, sigx_lb = 4, 6
sigtxr_lb = 8
r_lb = 10

ylims = [(-0.1, 0.1), (0.01, 300)]
img = Image_Plot()
img.set_style()
img.subplots(1, 1)

# img.axs[0][0].errorbar(result[r_lb], result[gt_lb], result[gt_lb + 1], c="C1", capsize=4, label="T", marker="s")
# img.axs[0][0].errorbar(result[r_lb], result[gx_lb], result[gx_lb + 1], c="C2", capsize=4, label="X", marker="s")

img.axs[0][0].errorbar(result[r_lb], result[sigt_lb], result[sigt_lb + 1], c="C1", capsize=4, label="T", marker="s")
img.axs[0][0].errorbar(result[r_lb], result[sigx_lb], result[sigx_lb + 1], c="C2", capsize=4, label="X", marker="s")


img.set_label(0, 0, 0, ylabels[1])
img.set_label(0, 0, 1, xlabel)

img.axs[0][0].set_yscale("log")
img.axs[0][0].set_ylim(ylims[1])

img.axs[0][0].set_xscale("log")
xs = img.axs[0][0].set_xlim()
img.axs[0][0].plot([xs[0], xs[1]], [0, 0], linestyle="--", linewidth=1, c="grey")
img.set_legend(0,0, loc="upper right")

for j in range(10):
    img.axs[0][0].plot([xs[0], xs[1]], [j, j], linewidth=0.7, c="grey", alpha=0.6)
    img.axs[0][0].plot([xs[0], xs[1]], [10 + 10*j, 10 + 10*j], linewidth=0.7,c="grey", alpha=0.6)
    img.axs[0][0].plot([xs[0], xs[1]], [100 + 100*j, 100 + 100*j], linewidth=0.7,c="grey", alpha=0.6)

img.axs[0][0].set_xlim(xs[0], xs[1])

img.save_img(dens_pic_path + ".png")
img.set_style_default()
img.close_img()

# plot R x \Delta\Sigma
img = Image_Plot()
img.set_style()
img.subplots(1,1)
img.set_label(0, 0, 0, ylabels_r)
img.set_label(0, 0, 1, xlabel)
img.axs[0][0].errorbar(result[r_lb], result[sigtxr_lb], result[sigtxr_lb + 1], c="C1", capsize=4, label="X", marker="s")
img.axs[0][0].set_xscale("log")
img.save_img(dens_r_pic_path + ".png")
img.set_style_default()
img.close_img()

print("Images are saved in %s"%dens_pic_path)


#
# for i in range(1):
#     img.set_label(0, i, 0, ylabels[i])
#     img.set_label(0, i, 1, xlabel)
#     if i == 1:
#         img.axs[0][i].set_yscale("log")
#         img.axs[0][i].set_ylim(ylims[i])
#
#     img.axs[0][i].set_xscale("log")
#     xs = img.axs[0][i].set_xlim()
#     img.axs[0][i].plot([xs[0], xs[1]], [0, 0], linestyle="--", linewidth=1, c="grey")
#     img.set_legend(0,i, loc="upper right")
#
#     if i == 1:
#         for j in range(10):
#             img.axs[0][i].plot([xs[0], xs[1]], [j, j], linewidth=0.7, c="grey", alpha=0.6)
#             img.axs[0][i].plot([xs[0], xs[1]], [10 + 10*j, 10 + 10*j], linewidth=0.7,c="grey", alpha=0.6)
#             img.axs[0][i].plot([xs[0], xs[1]], [100 + 100*j, 100 + 100*j], linewidth=0.7,c="grey", alpha=0.6)
#     img.axs[0][i].set_xlim(xs[0], xs[1])
# img.subimg_adjust(0,0.25)
# img.save_img(dens_pic_path + ".png")
# img.set_style_default()
# img.close_img()
#
# # plot R x \Delta\Sigma
# img = Image_Plot()
# img.set_style()
# img.subplots(1,1)
# img.set_label(0, 0, 0, ylabels_r)
# img.set_label(0, 0, 1, xlabel)
# img.axs[0][0].errorbar(result[r_lb], result[sigtxr_lb], result[sigtxr_lb + 1], c="C1", capsize=4, label="X", marker="s")
# img.axs[0][0].set_xscale("log")
# img.save_img(dens_r_pic_path + ".png")
# img.set_style_default()
# img.close_img()
#
# print("Images are saved in %s"%dens_pic_path)