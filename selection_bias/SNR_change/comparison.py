import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
from subprocess import Popen
import numpy
from plot_tool import Image_Plot




# plot
# matplotlib.style.use('default')
img = Image_Plot(fig_x=6, fig_y=4, ypad=0,xpad=0.2)
img.subplots(2, 1)
# img.set_style_default()
img.set_style()
img.axis_type(0, "major")
img.axis_type(1, "major")

for i in range(2):
    img.axs[i][0].set_ylabel("Relative error", fontsize=img.xy_lb_size + 2)

    img.axs[0][0].set_xticklabels([])

    img.axs[1][0].set_xlabel("$g_1$", fontsize=img.xy_lb_size + 2)


markers = ['o', 'v', 's', 'h', 'd', 'p', "4", "*", "X", "^", ">", "+"]
colors = ["C%d" % i for i in range(10)]



fmt = '%2.f%%'

flux_num = 4
parent_path = "/home/hklee/work/selection_bias/sensetive"

nms = ["/sex_pts/imgs_sex_sig/stack_result.npz",
       "/sex_gal/imgs_sex_sig/stack_result.npz"]

sig_label = ["Random walk galaxy\n\n$\\nu_F$",
             "Galsim galaxy\n\n$\\nu_F$"]

for i in range(2):
    m, n = divmod(i, 1)
    npz_file = numpy.load(parent_path + nms[i])
    shears = npz_file["arr_0"]
    snr_pool = npz_file["arr_5"]
    result_stack = npz_file["arr_1"][:flux_num]
    error_bar_stack = npz_file["arr_2"][: flux_num]

    for j in range(flux_num):
        idx_s = snr_pool[j] > 0

        snr_0 = snr_pool[j][idx_s].mean()
        var_rate = result_stack[j]
        var_err = error_bar_stack[j]
        if idx_s.sum() > 0:
            idx = var_rate > -5
            lb = "SNR$_t$ = %.2f"%snr_0
            # img.axs[m][n].scatter(shears[idx], var_rate[idx], edgecolor=colors[j], s=80, label=lb,
            #                       marker=markers[j], facecolor="none", linewidths=img.plt_line_width)
            img.axs[m][n].errorbar(shears[idx], var_rate[idx],var_err[idx], c=colors[j], ms=9,fmt=" ",label=lb,capsize=3,
                                   marker=markers[j], mfc="none", linewidth=img.plt_line_width)
            # print(labels[i],lb,var_rate)
x_ticks = numpy.linspace(-0.06, 0.06, 5)
for i in range(2):
    m,n = divmod(i,1)

    img.axs[m][n].set_ylim(-0.0035, 0.0035)
    img.axs[m][n].set_yticks(numpy.linspace(-0.003, 0.003, 5))

    img.axs[m][n].set_xlim(-0.075, 0.075)
    img.axs[m][n].set_xticks(x_ticks)

    img.axs_text(m, n, 0.8, 0.1,sig_label[i], text_color='k', text_fontsize=img.legend_size+1)
    img.axs[m][n].legend(fontsize=img.legend_size+1, loc="lower left", frameon=False,ncol=2,handletextpad=0.5)
# img.axs_text(m, n, text_y, text_x, labels[i], text_color='k', text_fontsize=img.legend_size+2)
img.save_img(parent_path + "/sigma_comparison.png")
img.save_img(parent_path + "/sigma_comparison.pdf")
img.close_img()
# img.show_img()