import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box


total_path = "E:/data/new_pdf/data_com"
h5f = h5py.File(total_path + "/shear.hdf5", "r")
g1 = h5f["/g1"][()]
g2 = h5f["/g2"][()]
h5f.close()
shear_tag = 3
print(g1[shear_tag], g2[shear_tag])
h5f = h5py.File(total_path + "/noise_free/data_%d_epsf.hdf5"%shear_tag, "r")
data = h5f["/data"][()]
mg1_noise_free = data[:,0]
mg2_noise_free = data[:,1]
mn_noise_free = data[:,2]
mu_noise_free = data[:,3]
h5f.close()


h5f = h5py.File(total_path + "/noisy/data_%d_epsf.hdf5"%shear_tag, "r")
data = h5f["/data"][()]
mg1_noisy = data[:,0]
mg2_noisy = data[:,1]
mn_noisy = data[:,2]
mu_noisy = data[:,3]
h5f.close()


h5f = h5py.File(total_path + "/noise_residual/data_%d_epsf.hdf5"%shear_tag, "r")
data = h5f["/data"][()]
mg1_noise_residual = data[:,0]
mg2_noise_residual = data[:,1]
mn_noise_residual = data[:,2]
mu_noise_residual = data[:,3]
h5f.close()

h5f = h5py.File(total_path + "/cross_term/data_%d_epsf.hdf5"%shear_tag, "r")
data = h5f["/data"][()]
mg1_cross_term = data[:,0]
mg2_cross_term = data[:,1]
mn_cross_term = data[:,2]
mu_cross_term = data[:,3]
h5f.close()

# d1 = mg1_noisy - mg1_noise_free - mg1_cross_term - mg1_noise_residual
# d2 = mg2_noisy - mg2_noise_free - mg2_cross_term - mg2_noise_residual
# d3 = mn_noisy - mn_noise_free - mn_cross_term - mn_noise_residual
# d4 = mu_noisy - mu_noise_free - mu_cross_term - mu_noise_residual
# print(d1.max(), d1.min())
# print(d2.max(), d2.min())
# print(d3.max(), d3.min())
# print(d4.max(), d4.min())
ls = 1.5
alpha=1
# # the hist of noise_free & noisy shear estimators
# img = Image_Plot(xpad=0.15, ypad=0.15)
# img.subplots(2,2)
# img.axs[0][0].hist(mg1_noise_free,100, label="noise_free G1",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][0].hist(mg1_noisy,100, label="noisy G1",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][0].hist(mg1_noise_residual,100, label="noise_residual G1",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][0].hist(mg1_cross_term,100, label="cross_term G1",alpha=alpha,histtype="step",linewidth=ls)
#
# img.axs[1][0].hist(mg2_noise_free,100, label="noise_free G2",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][0].hist(mg2_noisy,100, label="noisy G2",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][0].hist(mg2_noise_residual,100, label="noise_residual G2",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][0].hist(mg2_cross_term,100, label="cross_term G2",alpha=alpha,histtype="step",linewidth=ls)
#
# img.axs[0][1].hist(mn_noise_free,100, label="noise_free N",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][1].hist(mn_noisy,100, label="noisy N",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][1].hist(mn_noise_residual,100, label="noise_residual N",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][1].hist(mn_cross_term,100, label="cross_term N",alpha=alpha,histtype="step",linewidth=ls)
#
# img.axs[1][1].hist(mu_noise_free,100, label="noise_free U",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][1].hist(mu_noisy,100, label="noisy U",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][1].hist(mu_noise_residual,100, label="noise_residual U",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][1].hist(mu_cross_term,100, label="cross_term U",alpha=alpha,histtype="step",linewidth=ls)
#
# for i in range(2):
#     for j in range(2):
#         img.axs[i][j].legend()
#         img.axis_sci_ticklabel(i, j, 0)
#         img.axis_sci_ticklabel(i, j, 1)
# img.show_img()
#
# # ratio between the noisy & noise_free shear estimator
#
# img = Image_Plot(xpad=0.15, ypad=0.15)
# img.subplots(2,2)
# img.axs[0][0].hist(mg1_noise_free,100, label="noise_free G1",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][0].hist(mg1_noise_free + mg1_noise_residual,100, label="noise_free G1 + NR",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][0].hist(mg1_noise_free + mg1_cross_term,100, label="noise_free G1 + CT",alpha=alpha,histtype="step",linewidth=ls)
#
#
# img.axs[1][0].hist(mg2_noise_free,100, label="noise_free G2",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][0].hist(mg2_noise_free + mg2_noise_residual,100, label="noise_free G2 + NR",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][0].hist(mg2_noise_free + mg2_cross_term,100, label="noise_free G2 + CT",alpha=alpha,histtype="step",linewidth=ls)
#
#
# img.axs[0][1].hist(mn_noise_free,100, label="noise_free N",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][1].hist(mn_noise_free + mn_noise_residual,100, label="noise_free N + NR",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[0][1].hist(mn_noise_free + mn_cross_term,100, label="noise_free N + CT",alpha=alpha,histtype="step",linewidth=ls)
#
#
# img.axs[1][1].hist(mu_noise_free,100, label="noise_free U",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][1].hist(mu_noise_free + mu_noise_residual,100, label="noise_free U + NR",alpha=alpha,histtype="step",linewidth=ls)
# img.axs[1][1].hist(mu_noise_free + mu_cross_term,100, label="noise_free U + CT",alpha=alpha,histtype="step",linewidth=ls)
#
# for i in range(2):
#     for j in range(2):
#         img.axs[i][j].legend()
#         img.axis_sci_ticklabel(i, j, 0)
#         img.axis_sci_ticklabel(i, j, 1)
# img.show_img()

# idx1 = mg1_noise_free > 1.e3
# idx2 = mg1_noise_free < 1.e3 + 1000
# idx = idx1 & idx2
#
# plt_data = [[mg1_noise_free,mg1_noise_free + mg1_noise_residual,mg1_noise_free + mg1_cross_term],
#             [mn_noise_free, mn_noise_free + mn_noise_residual,mn_noise_free + mn_cross_term],
#             [mu_noise_free, mu_noise_free + mu_noise_residual,mu_noise_free + mu_cross_term],
#             [mn_noise_free+mu_noise_free, mn_noise_free + mn_noise_residual+mu_noise_free,
#              mn_noise_free + mn_cross_term+mu_noise_free]]
#
# labels = [["NF G1", "NF G1 + NR","NF G1 + CT"],["NF N", "NF N + NR", "NF N + CT"],
#           ["NF U","NF U + NR", "NF U + CT"], ["NF U + N", "NF U + N + NR", "NF U + N + CT"]]
#
# img = Image_Plot(fig_x=5, fig_y=4,xpad=0.1, ypad=0.13)
# img.subplots(2,2)
#
# for i in range(4):
#     m,n = divmod(i, 2)
#     img.axs[m][n].hist(plt_data[i][0][idx],20, label=labels[i][0],alpha=alpha,histtype="step",linewidth=ls)
#     img.axs[m][n].hist(plt_data[i][1][idx], 20, label=labels[i][1], alpha=alpha, histtype="step", linewidth=ls)
#     img.axs[m][n].hist(plt_data[i][2][idx],20, label=labels[i][2],alpha=alpha,histtype="step",linewidth=ls)
#
# for i in range(2):
#     for j in range(2):
#         img.axs[i][j].legend()
#         img.axis_sci_ticklabel(i, j, 0)
#         img.axis_sci_ticklabel(i, j, 1)
# img.show_img()


idx1 = mg2_noise_free > -1.e5
idx2 = mg2_noise_free < -1.e5 + 1000
idx = idx1 & idx2

plt_data = [[mg2_noise_free, mg2_noise_residual+mg2_noise_free, mg2_cross_term+mg2_noise_free],
            [mn_noise_free, mn_noise_residual, mn_cross_term],
            [mu_noise_free, mu_noise_residual, mu_cross_term],
            [mn_noise_free+mu_noise_free,  mn_noise_residual + mu_noise_free, mn_cross_term+mu_noise_free]]

labels = [["NF G1", "NR G1","CT G1"],["NF N", "NR N", "CT N"],
          ["NF U","NR U", "CT U"], ["NF U + N", "NR U + N", "CT U + N"]]

img = Image_Plot(fig_x=5, fig_y=4,xpad=0.1, ypad=0.13)
img.subplots(2,2)

for i in range(4):
    m,n = divmod(i, 2)
    img.axs[m][n].hist(plt_data[i][0][idx],20, label=labels[i][0],alpha=alpha,histtype="step",linewidth=ls)
    img.axs[m][n].hist(plt_data[i][1][idx], 20, label=labels[i][1], alpha=alpha, histtype="step", linewidth=ls)
    img.axs[m][n].hist(plt_data[i][2][idx],20, label=labels[i][2],alpha=alpha,histtype="step",linewidth=ls)

for i in range(2):
    for j in range(2):
        img.axs[i][j].legend()
        img.axis_sci_ticklabel(i, j, 0)
        img.axis_sci_ticklabel(i, j, 1)
img.show_img()
