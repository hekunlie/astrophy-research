import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
from subprocess import Popen
import numpy
from plot_tool import Image_Plot


cmd = argv[1]
# generate the galaxy & measurement
num = 300
if cmd == "run":
    for i in range(num):
        cmd = "python SNR_change.py 100 0.7 0.8 0.6 0.7 %d %d"%(1230+i, i)
        a = Popen(cmd, shell=True)
        a.wait()

# stack the result to give a mean measurement
if cmd == "stack":

    total_path = "./imgs/fits/"

    flux_num = 3
    snr_num = 4
    shear_num = 11

    # for the traditional SNR stacking
    snr_stack = numpy.zeros((int(flux_num*snr_num),))
    ori_detect = numpy.zeros((int(flux_num*snr_num),), dtype=numpy.intc)

    # for the result stacking
    result_stack = numpy.zeros((int(flux_num*snr_num), shear_num))
    detect_num = numpy.zeros((int(flux_num*snr_num), shear_num), dtype=numpy.intc)

    for k in range(num):
        # read the data
        file_npz = numpy.load(total_path + "result.npz")
        shears = file_npz["arr_0"]
        # only the 0 & 1'th column is needed
        ori_data, snr = file_npz["arr_1"][:, 0], file_npz["arr_1"][:, 1]
        shear_data = file_npz["arr_2"]

        for i in range(snr_num):
            for j in range(flux_num):
                # the non-detected is labeled by -99
                pre_shear_value, snr_0 = ori_data[j + i * flux_num], snr[j + i * flux_num]
                if pre_shear_value > 0:
                    sub_data = shear_data[j + i * flux_num]
                    # the detected
                    idx1 = sub_data > 0
                    # the non-detected
                    idx2 = sub_data < 0
                    var_rate = (sub_data - pre_shear_value) / pre_shear_value
                    var_rate[idx2] = 0

                    detection = numpy.ones((shear_num,), dtype=numpy.intc)
                    detection[idx2] = 0

                    result_stack[j + i*flux_num] += var_rate
                    detect_num[j + i*flux_num] += detection

                    # the original SNR
                    snr_stack[j + i*flux_num] += pre_shear_value
                    ori_detect[j + i*flux_num] += 1

    # calculate the mean
    idx = ori_detect < 1
    snr_stack[idx] = 0
    ori_detect[idx] = 1
    snr_stack = snr_stack/ori_detect

    idx = detect_num < 1
    detect_num[idx] = 1
    result_stack[idx] = 0







    result_stack = result_stack/detect_num

    numpy.savez("./imgs/stack_result.npz", shears,)

    # plot
    img = Image_Plot()
    img.subplots(1, snr_num)

    markers = ['o', 'v', 's', 'h', 'd', 'p', "4", "*", "X", "^", ">", "+"]
    colors = ["C%d" % i for i in range(10)]

    labels = ["P$_{k0}$", "SNR$_S$", "SNR$_A$", "MAG_AUTO"]

    fmt = '%2.f%%'

    for i in range(snr_num):
        if i == 0:
            img.axs[0][i].set_ylabel("Variation rate", fontsize=img.xy_lb_size)
        # else:
        #     img.axs[0][i].set_yticklabels([])

    for i in range(snr_num):
        for j in range(flux_num):
            # the non-detected is labeled by -99
            pre_shear_value, snr_0 = ori_data[j + i * flux_num], snr[j + i * flux_num]
            if pre_shear_value > 0:
                plt_data = shear_data[j + i * flux_num]
                idx = plt_data > 0
                var_rate = (plt_data - pre_shear_value) / pre_shear_value

                lb = "%s (%.2f)" % (labels[i], snr_0)
                img.axs[0][i].scatter(shears, var_rate, edgecolor=colors[j], s=80, label=lb,
                                      marker=markers[j], facecolor="none", linewidths=2)
            else:
                print("Non-detected. %s - %d" % (labels[i], j))
    # ys = [0,0]
    # for i in range(4):
    #     ys = img.axs[0][i].set_ylim()
    #     if ys[1] > ys[1]:
    #         ys[1] = ys[1]
    #     if ys[0] < ys[0]:
    #         ys[0] = ys[0]
    # dy = ys[1] - ys[0]
    # x_ticks = numpy.linspace(-0.06, 0.06, 5)
    # y_ticks = numpy.linspace(-0.04, 0.04, 5)

    for i in range(4):
        # img.axs[0][i].set_ylim(-0.042, 0.044)
        # img.axs[0][i].set_xlim(-0.075,0.075)
        # img.axs[0][i].set_xticks(x_ticks)
        # img.axs[0][i].set_yticks(y_ticks)
        img.axs[0][i].set_xlabel("$g_1$", fontsize=img.xy_lb_size)
        img.axs[0][i].legend(fontsize=img.legend_size, loc="best", frameon=False)

    # img.subimg_adjust(h=0, w=0)
    img.save_img(pic_path)
    print(pic_path)
    img.close_img()