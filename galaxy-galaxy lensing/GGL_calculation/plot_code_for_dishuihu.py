import numpy
import matplotlib.pyplot as plt
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
from Fourier_Quad import Fourier_Quad
import h5py


# cfht_parent_path = "F:/works/gg_lensing/cmass/cfht/6/"
# fourier_parent_path = "F:/works/gg_lensing/cmass/fourier/6/"
#
# ylabels = ["$\gamma$", "$\Delta\Sigma \; [\\rm{M_{\odot}} \cdot \\rm{pc^{-2}}]$"]
# ylabels_r = "$\\rm{R}\Delta\Sigma \; [\\rm{10^6\ M_{\odot}} \cdot \\rm{pc^{-2}}]$"
# xlabel = "$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$"
#
# fourier_names = ["fourier_cmass_result_total.hdf5", "fourier_cmass_result_w_1.hdf5", "fourier_cmass_result_w_3.hdf5"]
# cfht_names = ["cfht_cmass_result_total.hdf5", "cfht_cmass_result_w_1.hdf5", "cfht_cmass_result_w_3.hdf5"]
#
# # labels = ["corrected by <1+m>", "corrected by 1+\Sum(w_i m_i)/\Sum w_i"]
# labels = ["total", "$w_1$", "$w_3$"]
# colors=["C0", "C1", "C2"]
# lines = ["--", "-"]
#
# # CFHT cata
# img = Image_Plot()
# img.subplots(1,1)
# for i in range(3):
#     data_path = cfht_parent_path + "%s"%cfht_names[i]
#     h5f = h5py.File(data_path, "r")
#     result = h5f["/data"].value
#     h5f.close()
#     img.axs[0][0].errorbar(result[-1], result[4], result[5], linewidth=img.plt_line_width,capsize=3, c=colors[i],  label=labels[i])
#
# img.set_legend(0,0, loc="upper right")
# img.axs[0][0].set_xscale("log")
# img.axs[0][0].set_yscale("log")
#
# img.set_label(0, 0, 0, ylabels[1])
# img.set_label(0, 0, 1, xlabel)
# img.axs[0][0].set_ylim(0.1, 140)
# img.save_img("%sCFHT_cata.png"%cfht_parent_path)
# img.show_img()
# img.close_img()
#
# # Fourier_Quad cata
# img = Image_Plot()
# img.subplots(1, 1)
# for i in range(3):
#     data_path = fourier_parent_path + fourier_names[i]
#     h5f = h5py.File(data_path, "r")
#     result = h5f["/data"].value
#     h5f.close()
#     img.axs[0][0].errorbar(result[-1], result[4], result[5], linewidth=img.plt_line_width, capsize=3, c=colors[i], linestyle=lines[1], label=labels[i])
#
# img.set_legend(0,0, loc="upper right")
# img.axs[0][0].set_xscale("log")
# img.axs[0][0].set_yscale("log")
#
# img.set_label(0, 0, 0, ylabels[1])
# img.set_label(0, 0, 1, xlabel)
#
# img.save_img("%sfourier_cata.png"%fourier_parent_path)
# img.show_img()
# img.close_img()
#
#
# # comparison
# img = Image_Plot()
# img.subplots(1, 1)
# for i in range(1):
#
#     data_path = fourier_parent_path + fourier_names[i]
#     h5f = h5py.File(data_path, "r")
#     fourier_result = h5f["/data"].value
#     h5f.close()
#
#     data_path = cfht_parent_path + cfht_names[i]
#     h5f = h5py.File(data_path, "r")
#     cfht_result = h5f["/data"].value
#     h5f.close()
#
#     img.axs[0][0].errorbar(fourier_result[-1], fourier_result[4], fourier_result[5], capsize=3,
#                            linewidth=img.plt_line_width, c=colors[i], linestyle=lines[1], label="Fourier_Quad")
#
#     img.axs[0][0].errorbar(cfht_result[-1], cfht_result[4], cfht_result[5], capsize=3,
#                            linewidth=img.plt_line_width, c=colors[i+1], linestyle=lines[0], label=" CFHT")
# img.set_legend(0,0, loc="upper right")
# img.axs[0][0].set_xscale("log")
# img.axs[0][0].set_yscale("log")
# img.axs[0][0].set_ylim(0.1, 140)
# img.set_label(0, 0, 0, ylabels[1])
# img.set_label(0, 0, 1, xlabel)
# img.save_img("E:/works/dishuihu/compare.png")
# img.save_img("%sfourier_cata.png"%fourier_parent_path)
# img.show_img()
# img.close_img()
#
#
#
# # plot the comparison in each redshift bin of the Cluster catalog
# parent_path = "F:/works/gg_lensing/CFHT_cluster/"
#
# names = ["_CFHT_cluster_result_total.hdf5", "_CFHT_cluster_result_w_1.hdf5", "_CFHT_cluster_result_w_3.hdf5"]
#
# ylabels = ["$\gamma$", "$\Delta\Sigma \; [\\rm{M_{\odot}} \cdot \\rm{pc^{-2}}]$"]
# ylabels_r = "$\\rm{R}\Delta\Sigma \; [\\rm{10^6\ M_{\odot}} \cdot \\rm{pc^{-2}}]$"
# xlabel = "$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$"
#
# # 3 Z bins, each bin has 3 result: total, w1, w3
# cfht_data = [[], [], []]
# fourier_data = [[], [], []]
#
# z_bin = [0.2, 0.4, 0.6, 0.8]
# img = Image_Plot()
# img.subplots(1, 3)
# for i in range(3):
#     for j in range(3):
#         c_path = parent_path + 'cfht/6_%d/cfht%s'%(i+1, names[j])
#         h5f_c = h5py.File(c_path,"r")
#         data_c = h5f_c["/data"].value
#         h5f_c.close()
#         cfht_data[i].append(data_c)
#
#         f_path = parent_path + 'fourier/6_%d/fourier%s' % (i + 1, names[j])
#         h5f_f = h5py.File(f_path,"r")
#         data_f = h5f_f["/data"].value
#         h5f_f.close()
#         fourier_data[i].append(data_f)
#
#     plot_data_c = cfht_data[i][0]
#     plot_data_f = fourier_data[i][0]
#     img.axs[0][i].errorbar(plot_data_f[-1],plot_data_f[4],plot_data_f[5], capsize=3, c="C1", label="Fourier_Quad")
#     img.axs[0][i].errorbar(plot_data_c[-1],plot_data_c[4],plot_data_c[5], capsize=3, c="C2", label="CFHTLenS")
#     img.set_legend(0, i, loc="upper right")
#     img.axs[0][i].set_xscale("log")
#     img.axs[0][i].set_yscale("log")
#     img.set_label(0, i, 1, xlabel)
#     img.axs[0][i].set_ylim(0.1,140)
#     if i != 0:
#         img.axs[0][i].tick_params(labelleft=False)
#
#     text_cont = "$%.1f \leq Z < %.1f$"%(z_bin[i], z_bin[i+1])
#     img.axs[0][i].text(2, 20, text_cont, color='black', alpha=0.7, ha='left', va='center', fontsize=img.xy_lb_size)
#
# img.set_label(0, 0, 0, ylabels[1])
# img.subimg_adjust(0,0)
# img.save_img("E:/works/dishuihu/compare_.png")
# img.show_img()

# Overlap of CMASS and CFHT
# data_path = "F:/"
# area_range = numpy.array([[30.178211, 38.820651,-11.246086, -3.677158],
#                             [132.058306, 136.845897,-5.695229, -0.949215],
#                             [208.559162, 220.381114,51.196749, 57.805079],
#                             [329.977415, 335.710986,-1.029174, 4.614513]])
#
# ra_min = area_range[:,0]
# ra_max = area_range[:,1]
# dec_min = area_range[:,2]
# dec_max = area_range[:,3]
#
# img = Image_Plot()
# img.subplots(2,2)
# img.set_style()
# redshift = []
# for i in range(4):
#     cmass_path = data_path + "w_%d.hdf5"%(i+1)
#     h5f = h5py.File(cmass_path,"r")
#     ra = h5f["/RA"].value
#     dec = h5f["/DEC"].value
#     z = h5f["/Z"].value
#     if i != 1:
#         redshift.extend(z.tolist())
#     h5f.close()
#     m,n = divmod(i,2)
#
#     norm = plt.Normalize(vmin=numpy.min(z), vmax=numpy.max(z))
#     cmap = plt.get_cmap('jet')
#     cl = cmap(norm(z))
#
#     idx_1 = ra <= ra_max[i]
#     idx_2 = ra >= ra_min[i]
#     idx_3 = dec <= dec_max[i]
#     idx_4 = dec >= dec_min[i]
#
#     idx = idx_1 & idx_2 & idx_3 & idx_4
#
#     img.axs[m][n].scatter(ra[idx], dec[idx],color=cl)
#     sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
#     sm._A = []
#
#     plt.colorbar(sm, ax=img.axs[m][n])
#     # img.axs[m][n].plot([ra_min[0], ra_min[0]],[ra_min[0], ra_min[0]])
#     img.set_label(m, n, 0, "DEC. [deg]")
#     img.set_label(m, n, 1, "RA. [deg]")
# img.subimg_adjust(0.2, 0.2)
# img.save_img("E:/works/dishuihu/cmass.png")
# img.show_img()
# img.close_img()
#
# print(len(redshift))
# img = Image_Plot()
# img.subplots(1,1)
# img.set_style()
# img.axs[0][0].hist(redshift,20,color="blue")
# img.set_label(0,0, 1, "$Z$")
# img.set_label(0, 0, 0, "$N(Z)$")
# img.save_img("E:/works/dishuihu/cmass_PDF_Z.png")
# img.show_img()
# img.close_img()

# # cluster N200
# data_path = "E:/works/dishuihu/cluster/"
#
# z_bin = [0.2, 0.4, 0.6, 0.8]
# img = Image_Plot()
# img.subplots(1,3)
# img.set_style()
# ys = []
# for i in range(3):
#     for j in range(4):
#         h5f = h5py.File(data_path + "w_%d.hdf5"%(j+1), "r")
#         temp = h5f["/Z%d/N200"%i].value
#         temp.shape = (temp.shape[0],1)
#         h5f.close()
#         if j == 0:
#             n200 = temp
#         else:
#             n200 = numpy.row_stack((n200, temp))
#     img.axs[0][i].hist(n200[n200<100],40,color="blue")
#
#     img.set_label(0, i, 1, "$N200$")
#
#     text_cont = "$%.1f \leq Z < %.1f$" % (z_bin[i], z_bin[i + 1])
#     img.axs[0][i].text(50, 800, text_cont, color='black', alpha=0.7, ha='left', va='center', fontsize=img.xy_lb_size)
#
#     if i != 0:
#         img.axs[0][i].tick_params(labelleft=False)
#     ys_ = img.axs[0][i].set_ylim()
#     ys.append(ys_[1])
#
# for i in range(3):
#     img.axs[0][i].set_ylim(0,max(ys))
# img.set_label(0, 0, 0, "$N$")
# img.subimg_adjust(0,0)
# img.save_img("E:/works/dishuihu/cluster_PDF_N200.png")
# img.show_img()
# img.close_img()

size = 48
sigma = 4
sig_level = 3

cen = int(size/2)
nx, ny = size, size
ra_bin = numpy.linspace(0,size, size+1)
dec_bin = ra_bin
ra_min, ra_max = ra_bin.min(), ra_bin.max()
dec_min, dec_max = dec_bin.min(), dec_bin.max()

kappa = tool_box.gauss_profile(size, sigma, cen, cen)
shear = tool_box.kappa2shear(kappa)

gamma1 = shear.real
gamma2 = shear.imag

g1_noise = numpy.random.normal(0, numpy.abs(gamma1)/sig_level)
g2_noise = numpy.random.normal(0, numpy.abs(gamma2)/sig_level)

gamma1 = gamma1 + g1_noise
gamma2 = gamma2 + g2_noise

kappa_recon = tool_box.shear2kappa(gamma1, gamma2)
kappa_new = kappa_recon.real
kappa_curl = kappa_recon.imag

kappa_recon_n = tool_box.shear2kappa(gamma1, gamma2)
kappa_new_n = kappa_recon_n.real
kappa_curl_n = kappa_recon_n.imag

g = numpy.abs(shear)
g = numpy.sqrt(gamma1**2 + gamma2**2)
npw = numpy.where(g == 0)
cos_theta = numpy.sqrt((1+gamma1/g)/2)
sin_theta = gamma2/2/g/cos_theta
idx = cos_theta == 0
sin_theta[idx] = 1
max_g = g.max()

ra_sq_len = ra_bin[2] - ra_bin[1]
dec_sq_len = ra_sq_len
max_len = ra_sq_len *1.2



# img = Image_Plot()
# img.subplots(1,1)
# img.set_style()
#
# fig = img.axs[0][0].imshow(kappa)
# img.axs[0][0].set_title("$\kappa$ field", fontsize=img.xy_lb_size)
# plt.colorbar(fig, ax=img.axs[0][0])
# img.axs[0][0].tick_params(labelleft=False)
# img.axs[0][0].tick_params(labelbottom=False)
#
# img.save_img("E:/works/dishuihu/kappa_field.png")
# img.show_img()
# img.close_img()

img = Image_Plot(fig_x=8,fig_y=8)
img.subplots(1,1)
img.set_style()

img.axs[0][0].set_title("$\gamma$ field", fontsize=img.xy_lb_size)

for i in range(ny + 1):
    img.axs[0][0].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--", alpha=0.6, linewidth=0.3)
for j in range(nx + 1):
    img.axs[0][0].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black", linestyle="--", alpha=0.6, linewidth=0.3)


dg_scale = g / max_g * max_len / 2


img.axs[0][0].tick_params(labelleft=False)
img.axs[0][0].tick_params(labelbottom=False)
for i in range(ny):
    for j in range(nx):

        if gamma2[i, j] < 0:
            dx = -dg_scale[i, j] * cos_theta[i, j]
            dy = -dg_scale[i, j] * sin_theta[i, j]

        else:
            dx = dg_scale[i, j] * cos_theta[i, j]
            dy = dg_scale[i, j] * sin_theta[i, j]

        x = (ra_bin[j] + ra_bin[j + 1]) / 2
        y = (dec_bin[i] + dec_bin[i + 1]) / 2
        img.axs[0][0].plot([x + dx, x - dx], [y + dy, y - dy], c="C1")

img.save_img("E:/works/dishuihu/shear_field.png")
img.show_img()
img.close_img()




img = Image_Plot()
img.subplots(1, 1)
img.set_style()

fig = img.axs[0][0].imshow(kappa_new)
img.axs[0][0].set_title("$\kappa$ field", fontsize=img.xy_lb_size)
plt.colorbar(fig, ax=img.axs[0][0])
img.axs[0][0].tick_params(labelleft=False)
img.axs[0][0].tick_params(labelbottom=False)
img.save_img("E:/works/dishuihu/shear_kappa_recovery.png")
img.show_img()
img.close_img()

img = Image_Plot()
img.subplots(1, 1)
img.set_style()

fig = img.axs[0][0].imshow(kappa_new_n)
img.axs[0][0].set_title("$\kappa$ field", fontsize=img.xy_lb_size)
plt.colorbar(fig, ax=img.axs[0][0])
img.axs[0][0].tick_params(labelleft=False)
img.axs[0][0].tick_params(labelbottom=False)
img.save_img("E:/works/dishuihu/shear_kappa_recovery_n.png")
img.show_img()
img.close_img()
