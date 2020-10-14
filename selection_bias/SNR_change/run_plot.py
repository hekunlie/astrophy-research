import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
from subprocess import Popen
import numpy
from plot_tool import Image_Plot


source_label = argv[1]
cmd = argv[2]
# generate the galaxy & measurement
num_st = int(argv[3])
num_ed = int(argv[4])

if cmd == "run":
    for i in range(num_st, num_ed):
        cmd = "python SNR_change.py 100 0.7 0.7 0.6 0.7 %d %d %s"%(122230+i, i, source_label)
        a = Popen(cmd, shell=True)
        a.wait()

# stack the result to give a mean measurement
if cmd == "stack":

    total_path = "./imgs/fits/"
    pic_path_png = "./imgs/stack_result.png"
    pic_path_pdf = "./imgs/stack_result.pdf"

    flux_num = 4
    snr_num = 4
    shear_num = 11

    # for the traditional SNR stacking
    snr_stack = numpy.zeros((int(flux_num*snr_num), ))
    ori_detect = numpy.zeros((int(flux_num*snr_num), ), dtype=numpy.intc)

    # pre-sheared galaxy detection
    pre_shear_detect = numpy.zeros((int(flux_num*snr_num), ), dtype=numpy.intc)

    # for the result stacking
    result_stack = numpy.zeros((int(flux_num*snr_num), shear_num))
    error_bar_stack = numpy.zeros((int(flux_num*snr_num), shear_num))
    detect_num = numpy.zeros((int(flux_num*snr_num), shear_num), dtype=numpy.intc)

    # append each result to the pool for plot
    pk_pool = numpy.zeros((flux_num, num_ed-num_st,shear_num)) - 999
    snr_sex_pool = numpy.zeros((flux_num, num_ed-num_st,shear_num)) - 999
    mag_auto_pool = numpy.zeros((flux_num, num_ed-num_st,shear_num)) - 999
    rfacotr_pool = numpy.zeros((flux_num, num_ed-num_st,shear_num)) - 999
    snr_pool = numpy.zeros((flux_num, num_ed-num_st)) - 999

    for k in range(num_st,num_ed):
        # read the data
        file_npz = numpy.load(total_path + "%d/result.npz"%k)
        shears = file_npz["arr_0"]
        # only the 0 & 1'th column is needed
        ori_data, snr = file_npz["arr_1"][:, 0], file_npz["arr_1"][:, 1]
        # the "SNR's" after shear
        shear_data = file_npz["arr_2"]

        for i in range(flux_num):
            ori_pk0, ori_snr = ori_data[i], snr[i]
            pk0 = shear_data[i]

            ori_snr_sex = ori_data[i + flux_num]
            snr_sex = shear_data[i + flux_num]

            ori_mag_auto = ori_data[i + int(3*flux_num)]
            mag_auto = shear_data[i + int(3*flux_num)]

            ori_rfactor = ori_data[i + int(4*flux_num)]
            rfactor = shear_data[i + int(4*flux_num)]
            # print(ori_snr)
            if ori_snr > 0:
                snr_pool[i,k-num_st] = ori_snr

                idx = pk0 > 0
                pk_pool[i][k-num_st][idx] = (pk0[idx] - ori_pk0)/ori_pk0

                idx = snr_sex > 0
                snr_sex_pool[i][k-num_st][idx] = (snr_sex[idx] - ori_snr_sex)/ori_snr_sex

                idx = mag_auto > 0
                mag_auto_pool[i][k-num_st][idx] = (mag_auto[idx] - ori_mag_auto)/ori_mag_auto

                idx = rfactor > 0
                rfacotr_pool[i][k-num_st][idx] = (rfactor[idx] - ori_rfactor)/ori_rfactor

    for i in range(flux_num):
        for j in range(shear_num):
            idx = pk_pool[i][:,j] > -5
            result_stack[i, j] = pk_pool[i][:,j][idx].mean()
            error_bar_stack[i, j] = pk_pool[i][:,j][idx].std()

            idx = mag_auto_pool[i][:,j] > -5
            result_stack[i + flux_num, j] = mag_auto_pool[i][:,j][idx].mean()
            error_bar_stack[i + flux_num, j] = mag_auto_pool[i][:,j][idx].std()

            idx = snr_sex_pool[i][:,j] > -5
            result_stack[i + int(2*flux_num), j] = snr_sex_pool[i][:,j][idx].mean()
            error_bar_stack[i + int(2*flux_num), j] = snr_sex_pool[i][:,j][idx].std()

            idx = rfacotr_pool[i][:,j] > -5
            result_stack[i + int(3*flux_num), j] = rfacotr_pool[i][:,j][idx].mean()
            error_bar_stack[i + int(3*flux_num), j] = rfacotr_pool[i][:,j][idx].std()

    print(result_stack)
    #     for i in range(snr_num):
    #         for j in range(flux_num):
    #             # the non-detected is labeled by -99
    #             pre_shear_value, snr_0 = ori_data[j + i * flux_num], snr[j + i * flux_num]
    #             if pre_shear_value > 0:
    #                 # print(i,j,"Found",pre_shear_value)
    #                 sub_data = shear_data[j + i * flux_num]
    #                 # the non-detected
    #                 idx2 = sub_data <= 0
    #                 var_rate = (sub_data - pre_shear_value) / pre_shear_value
    #                 var_rate[idx2] = 0
    #
    #                 detection = numpy.ones((shear_num,), dtype=numpy.intc)
    #                 detection[idx2] = 0
    #                 # print(detection)
    #
    #                 result_stack[j + i * flux_num] += var_rate
    #                 detect_num[j + i * flux_num] += detection
    #
    #                 # the original SNR
    #                 snr_stack[j + i * flux_num] += snr_0
    #                 ori_detect[j + i * flux_num] += 1
    #                 # print(k, i, j, snr_0)
    #     print(k, ori_detect)
    #
    # # calculate the mean
    # print(snr_stack)
    # print(ori_detect)
    # idx = ori_detect < 1
    # ori_detect[idx] = 1
    # snr_stack = snr_stack/ori_detect
    # snr_stack[idx] = -1000
    #
    # idx = detect_num < 1
    # detect_num[idx] = 1
    # result_stack = result_stack/detect_num
    # result_stack[idx] = -1000
    # print(result_stack)


    numpy.savez("./imgs/stack_result.npz", shears, result_stack,error_bar_stack, snr_stack)

    # plot
    img = Image_Plot(fig_x=6, fig_y=4, ypad=0.22,xpad=0)
    img.subplots(2, 2)
    img.set_style_default()
    img.set_style()
    img.axis_type(0, "major")
    img.axis_type(1, "major")
    markers = ['o', 'v', 's', 'h', 'd', 'p', "4", "*", "X", "^", ">", "+"]
    colors = ["C%d" % i for i in range(10)]
    plot_data = [result_stack[: flux_num],
                 result_stack[flux_num: 2*flux_num],
                 result_stack[2*flux_num: 3*flux_num],
                 result_stack[3*flux_num: 4*flux_num]]
    plot_data_err = [error_bar_stack[: flux_num],
                    error_bar_stack[flux_num: 2*flux_num],
                    error_bar_stack[2*flux_num: 3*flux_num],
                    error_bar_stack[3*flux_num: 4*flux_num]]
    labels = ["P$_{k0}}$",  "MAG_AUTO", "SNR$_S$", "Resolution"]

    fmt = '%2.f%%'

    for i in range(len(labels)):

        m, n = divmod(i, 2)

        if i==0 or i==2:
            img.axs[m][n].set_ylabel("Relative error", fontsize=img.xy_lb_size+2)
        else:
            img.axs[m][n].set_yticklabels([])

        img.axs[m][n].set_xlabel("$g_1$", fontsize=img.xy_lb_size+2)


    for i in range(len(labels)):

        m, n = divmod(i, 2)
        for j in range(flux_num):
            idx_s = snr_pool[j] > 0

            snr_0 = snr_pool[j][idx_s].mean()
            var_rate = plot_data[i][j]
            var_err = plot_data_err[i][j]
            if idx_s.sum() > 0:
                idx = var_rate > -5
                lb = "SNR_i = %.2f"%snr_0
                # img.axs[m][n].scatter(shears[idx], var_rate[idx], edgecolor=colors[j], s=80, label=lb,
                #                       marker=markers[j], facecolor="none", linewidths=img.plt_line_width)
                img.axs[m][n].errorbar(shears[idx], var_rate[idx],var_err[idx], c=colors[j], ms=6,fmt=" ",label=lb,capsize=3,
                                       marker=markers[j], mfc="none", linewidth=img.plt_line_width)
                print(var_rate)
                print(var_err)
        if i == 0:
            text_x, text_y = 0.8, 0.1
        elif i == 1:
            text_x, text_y = 0.68, 0.1
        elif i == 2:
            text_x, text_y = 0.8, 0.9
        else:
            text_x, text_y = 0.68, 0.9
        img.axs_text(m, n, text_y, text_x, labels[i], text_color='k', text_fontsize=img.legend_size+2)
    ys = [0,0]
    for i in range(4):
        m, n = divmod(i, 2)
        ys_ = img.axs[m][n].set_ylim()
        if ys_[1] > ys[1]:
            ys[1] = ys_[1]
        if ys_[0] < ys[0]:
            ys[0] = ys_[0]
    x_ticks = numpy.linspace(-0.06, 0.06, 5)

    for i in range(4):
        m,n = divmod(i,2)
        if i < 2:
            if source_label == "galsim":
                img.axs[m][n].set_ylim(-0.0025, 0.0025)
                img.axs[m][n].set_yticks(numpy.linspace(-0.002, 0.002, 5))
            else:
                img.axs[m][n].set_ylim(-0.0025, 0.0025)
                img.axs[m][n].set_yticks(numpy.linspace(-0.002, 0.002, 5))
        else:
            img.axs[m][n].set_ylim(-0.037, 0.037)
            img.axs[m][n].set_yticks(numpy.linspace(-0.03, 0.03, 5))
        img.axs[m][n].set_xlim(-0.075, 0.075)
        img.axs[m][n].set_xticks(x_ticks)

        if i == 0:
            img.axs[m][n].legend(fontsize=img.legend_size+1, loc="upper left", frameon=False,ncol=2)

    # img.subimg_adjust(h=0, w=0)
    img.save_img(pic_path_png)
    img.save_img(pic_path_pdf)
    print(pic_path_pdf)
    img.close_img()