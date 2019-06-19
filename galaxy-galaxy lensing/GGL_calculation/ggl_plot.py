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

total_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/"

h5f = h5py.File(total_path + "data/cata_result_ext_grid.hdf5", "r")
radius_bin = h5f["/radius_bin"].value
sp = radius_bin.shape

radius = (radius_bin[1:] + radius_bin[:sp[0]-1])/2
h5f.close()

h5f_g = h5py.File(total_path + "result/%s/%s/result_gamma.hdf5"%(fore_source, data_file), "r")
gamma = h5f_g["/result"].value*100
chi_g = h5f_g["/chisq"].value # (13, 40)
h5f_g.close()

h5f_c = h5py.File(total_path + "result/%s/%s/result_crit.hdf5"%(fore_source, data_file), "r")
crit = h5f_c["/result"].value
chi_c = h5f_c["/chisq"].value
h5f_c.close()

signals = [gamma, crit]
chisq = [chi_g, chi_c]
ylims = [(-1.5,3.1), (-20, 300)]
ylabels = ["$\gamma \\times 10^2$", "$\Delta\Sigma$"]

img = Image_Plot()
img.subplots(1,2)
img.set_style()
for i in range(2):
    img.axs[0][i].errorbar(radius, signals[i][:, 0], signals[i][:, 1], c="C1", capsize=4, label="T", marker="p")
    img.axs[0][i].errorbar(radius, signals[i][:, 2], signals[i][:, 3], c="C2", capsize=4, label="X", marker="p")

    img.set_label(0,i,0, ylabels[i])
    img.set_label(0,i,1, "$R \quad \\rm{Mpc \cdot h^{-1}}$")
    if i == 1:
        img.axs[0][i].set_yscale("symlog")
    img.axs[0][i].set_xscale("log")

    img.axs[0][i].set_ylim(ylims[i])
    xs = img.axs[0][i].set_xlim()
    img.axs[0][i].plot([xs[0], xs[1]], [0, 0], linestyle="--", linewidth=1, c="grey")
    img.set_legend(0,i)

    if i == 1:
        for j in range(15):
            if j < 10:
                img.axs[0][i].plot([xs[0], xs[1]], [j, j], linewidth=1, c="grey", alpha=0.6)
            img.axs[0][i].plot([xs[0], xs[1]], [10 + 10*j, 10 + 10*j], linewidth=1,c="grey", alpha=0.6)
    img.axs[0][i].set_xlim(xs[0], xs[1])

img.save_img(total_path + "result/%s/%s/result_%s.pdf"%(fore_source, data_file, data_file))
img.save_img(total_path + "result/%s/%s/result_%s.png"%(fore_source, data_file, data_file))
img.close_img()

# img = Image_Plot(fig_y=3, fig_x=4)
# img.subplots(2,13)
# img.set_style()
# for i in range(2):
#     for j in range(13):
#         chi1 = chisq[i][j,:20]
#         chi2 = chisq[i][j,20:]
#         img.axs[i][j].scatter(chisq[i])

