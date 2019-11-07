import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import plot_tool
from mpi4py import MPI
import tool_box
import numpy
from sys import argv
import h5py
from plot_tool import Image_Plot
import matplotlib.ticker as mtick


# This program plots the figure of each cutoff

def mcplot(x, y, y_err, y_num, mc, cut_start, xylim, fig_ax):
    # "x_data' is the 'x'
    # 'y_data' is an (3,n) array "[[y's],[dy's],[num's]]

    fig_ax.errorbar(x, y, y_err, ecolor='C1', fmt=' ', capsize=4)
    fig_ax.plot(x, (mc[0]+1) * x + mc[2], color='k')
    fig_ax.plot([-0.1, 0.1], [-0.1, 0.1], label='y=x', color='blue')
    fig_ax.scatter(x, y, c='black')

    for j in range(len(x)):
        fig_ax.text(x[j], y[j], str(round(y_num[j] / 1000, 1)) + "K", color="red")
    fig_ax.text(0.1, 0.85, 'm=' + str(round(mc[0], 6)) + '$\pm$' + str(round(mc[1], 6)), color='green', ha='left',
            va='center', transform=fig_ax.transAxes, fontsize=20)
    fig_ax.text(0.1, 0.8, 'c=' + str(round(mc[2], 6)) + '$\pm$' + str(round(mc[3], 6)), color='green', ha='left',
            va='center', transform=fig_ax.transAxes, fontsize=20)
    fig_ax.text(0.1, 0.75, " $\geq$ %.5f"%cut_start, color='green', ha='left', va='center',
                transform=fig_ax.transAxes, fontsize=20)
    fig_ax.set_xlabel('True  g', fontsize=20)
    fig_ax.set_ylabel('Est  g', fontsize=20)
    fig_ax.legend(fontsize=15)
    fig_ax.set_ylim(xylim[0], xylim[1])
    fig_ax.set_xlim(xylim[0], xylim[1])


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

# the number of thread must be equal to cutoff number
total_path = argv[1]
filter_name = argv[2]
select_name = argv[3]


shear_path = total_path + "/parameters/shear.npz"
g1 = numpy.load(shear_path)["arr_0"]
g2 = numpy.load(shear_path)["arr_1"]

# range of plot
g1_extr = [g1.min(), g1.max()]
g2_extr = [g2.min(), g2.max()]
g1_plt_s, g1_plt_e = g1.min() - (g1.max() - g1.min()) * 0.2, g1.max() + (g1.max() - g1.min()) * 0.2
g2_plt_s, g2_plt_e = g2.min() - (g2.max() - g2.min()) * 0.2, g2.max() + (g2.max() - g2.min()) * 0.2
plt_range = [g1_plt_s, g1_plt_e, g2_plt_s, g2_plt_e]


file_path = total_path + "/result/cuts/sym/%s/%s/"%(filter_name, select_name)
h5f = h5py.File(file_path + "total.hdf5","r")

mc1 = h5f["/mc1"].value
mc2 = h5f["/mc2"].value
results = h5f["/shear"].value
num = h5f["/num"].value
cut_scale = h5f["/cut_scale"].value
h5f.close()

# the row labels are corresponding to the shears
shear_num, cut_num = mc1.shape


mc_title = ['0', '0', '0', '0']

e1mc = (mc1[0, rank], mc1[1, rank], mc1[2, rank], mc1[3, rank])
e2mc = (mc2[0, rank], mc2[1, rank], mc2[2, rank], mc2[3, rank])

# m1/2 +- 2sig
m_r = [[e1mc[0] - 2 * e1mc[1], e1mc[0] + 2 * e1mc[1]],
       [e2mc[0] - 2 * e2mc[1], e2mc[0] + 2 * e2mc[1]]]

c_r = [[e1mc[2] - 2 * e1mc[3], e1mc[2] + 2 * e1mc[3]],
       [e2mc[2] - 2 * e2mc[3], e2mc[2] + 2 * e2mc[3]]]

for ii in range(2):
    if tool_box.check_in(m_r[ii]):
        mc_title[ii] = ''
    else:
        mc_title[ii] = "_m" + str(ii+1)
    if tool_box.check_in(c_r[ii]):
        mc_title[ii + 2] = ''
    else:
        mc_title[ii + 2] = "_c" + str(ii+1)
pic_mc = "".join(mc_title)

pic_path = file_path + "cut_%d%s.png"%(rank, pic_mc)

# if rank == 0:
#     print(results)
#     print(results.shape)
img = Image_Plot(fig_x=10,fig_y=10)
img.subplots(1,2)

mcplot(g1, results[:,rank], results[:,rank+cut_num], num[:,rank], e1mc,
       cut_scale[rank],  plt_range[:2], img.axs[0][0])

mcplot(g2, results[:,rank+2*cut_num], results[:,rank+3*cut_num], num[:,rank], e2mc,
       cut_scale[rank], plt_range[2:4], img.axs[0][1])

img.save_img(pic_path)
img.close_img()

if rank == 0:

    ylabels = ["$m$", "$c$"]

    pic_path = file_path + "total.png"
    ch_num = 10
    x_coord = [i * cut_num for i in range(ch_num)]
    ch = [i for i in range(ch_num)]

    fmt = '%2.f%%'
    xticks = mtick.FormatStrFormatter(fmt)
    img = Image_Plot(fig_x=10,fig_y=10)
    img.subplots(1,2)

    img.axs[0][0].errorbar(x_coord, mc1[0], mc1[1], c="C1", linewidth=img.plt_line_width,capsize=img.cap_size, label="$m_1$", marker="s")
    img.axs[0][0].errorbar(x_coord, mc2[0], mc2[1], c="C2", linewidth=img.plt_line_width, capsize=img.cap_size, label="$m_2$", marker="s")

    img.axs[0][1].errorbar(x_coord, mc1[2], mc1[3], c="C1", linewidth=img.plt_line_width, capsize=img.cap_size, label="$c_1$", marker="s")
    img.axs[0][1].errorbar(x_coord, mc2[2], mc2[3], c="C2", linewidth=img.plt_line_width, capsize=img.cap_size, label="$c_2$", marker="s")

    for i in range(2):
        xs = img.axs[0][i].set_xlim()
        img.axs[0][i].plot([xs[0], 100], [0, 0], linewidth=img.plt_line_width, c="grey", linestyle="--")
        img.axs[0][i].set_xlim(xs[0], xs[1])
        img.axs[0][i].xaxis.set_major_formatter(xticks)
        img.axs[0][i].set_xlabel("Cutoff percentage", fontsize=img.xy_lb_size)
        img.axs[0][i].set_ylabel(ylabels[i], fontsize=img.xy_lb_size)
        img.axs[0][i].legend(fontsize=img.legend_size)
    img.save_img(pic_path)

# filter_name = ["sex2_4", "sex2_2", "sex2_1.5", "sex4_4", "sex4_2", "sex4_1.5"]
# select_name = ["mag_auto", "snr_auto", 'snr_sex', "flux2_ex1", "flux2_ex2", "flux2_ex3", "rfactor"]

# file_paths = []
# for fnm in filter_name:
#     for snm in select_name:
#         file_path = total_path + "result/cuts/sym/%s/%s/"%(fnm, snm)
#         file_paths.append(file_path)
#
# my_files = tool_box.allot(file_paths, cpus)[rank]
#
# for file_path in my_files:
#
#     pdfs = os.listdir(file_path)
#     for pdf_file in pdfs:
#         if ".pdf" in pdf_file:
#             os.remove(file_path + pdf_file)
#     h5f = h5py.File(file_path + "total.hdf5","r")
#
#     mc1 = h5f["/mc1"].value
#     mc2 = h5f["/mc2"].value
#     results = h5f["/shear"].value
#     num = h5f["/num"].value
#     cut_scale = h5f["/cut_scale"].value
#
#     # the row labels are corresponding to the shears
#     shear_num, cut_num = mc1.shape
#
#     h5f.close()
#
#     for i in range(cut_num):
#         mc_title = ['0', '0', '0', '0']
#         e1mc = (mc1[0, i]+1, mc1[1, i], mc1[2, i], mc1[3, i])
#         e2mc = (mc2[0, i]+1, mc2[1, i], mc2[2, i], mc2[3, i])
#         # m1/2 +- 2sig
#         m_r = [[e1mc[0] - 1 - 2 * e1mc[1], e1mc[0] - 1 + 2 * e1mc[1]],
#                [e2mc[0] - 1 - 2 * e2mc[1], e2mc[0] - 1 + 2 * e2mc[1]]]
#         c_r = [[e1mc[2] - 2 * e1mc[3], e1mc[2] + 2 * e1mc[3]],
#                [e2mc[2] - 2 * e2mc[3], e2mc[2] + 2 * e2mc[3]]]
#
#         for ii in range(2):
#             if tool_box.check_in(m_r[ii]):
#                 mc_title[ii] = ''
#             else:
#                 mc_title[ii] = "_m" + str(ii+1)
#             if tool_box.check_in(c_r[ii]):
#                 mc_title[ii + 2] = ''
#             else:
#                 mc_title[ii + 2] = "_c" + str(ii+1)
#         pic_mc = "".join(mc_title)
#
#         pic_path = file_path + "/" + str(round(cut_scale[0,i], 4)) + pic_mc+".pdf"
#
#         tool_box.mcplot(g1, results[:,i], results[:,i+cut_num], num[:,i],
#                         g2, results[:,i+2*cut_num], results[:,i+3*cut_num], num[:,i],
#                         e1mc, e2mc, str(round(cut_scale[0,i],4)), 'max', plt_range,pic_path)




