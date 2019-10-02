import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
import matplotlib.pyplot as plt
import tool_box
from plot_tool import Image_Plot


# subplot num, flux num, the data path and the image path
img_num, flux_num, data_path, pic_path = int(argv[1]), int(argv[2]), argv[3], argv[4]

file_tag = 1
matplotlib.style.use('default')
plt.rcParams['font.family'] = 'serif'

total_path = "./imgs/"

# read the data
file_npz = numpy.load(data_path)
shears = file_npz["arr_0"]
# only the 0 & 1'th column is needed
ori_data, snr = file_npz["arr_1"][:,0],file_npz["arr_1"][:,1]
shear_data = file_npz["arr_2"]

# plot
img = Image_Plot()
img.subplots(1, img_num)

markers = ['o','v','s','h','d','p',"4","*","X","^",">","+"]
colors = ["C%d"%i for i in range(10)]

labels = ["P$_{k0}$", "SNR$_S$", "SNR$_A$", "MAG_AUTO"]

fmt = '%2.f%%'

for i in range(img_num):
    if i == 0:
        img.axs[0][i].set_ylabel("Variation rate",fontsize=img.xy_lb_size)
    # else:
    #     img.axs[0][i].set_yticklabels([])


for i in range(img_num):
    for j in range(flux_num):
        # the non-detected is labeled by -99
        pre_shear_value, snr_0 = ori_data[j + i*flux_num], snr[j+i*flux_num]
        if pre_shear_value > 0:
            plt_data = shear_data[j + i*flux_num]
            idx = plt_data > 0
            var_rate = (plt_data - pre_shear_value)/pre_shear_value

            lb = "%s (%.2f)"%(labels[i], snr_0)
            # axs[select].plot(numpy.linspace(-0.06, 0.06, num)[idx], deltas[select][idx], c=colors[i], ms=12, label=lb,
            #          marker=markers[i],linestyle=' ',fillstyle='none')
            img.axs[0][i].scatter(shears[idx], var_rate[idx], edgecolor=colors[j], s=80, label=lb,
                     marker=markers[j], facecolor="none", linewidths=2)
        else:
            print("Non-detected. %s - %d"%(labels[i],j))
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
    img.axs[0][i].set_xlabel("$g_1$",fontsize=img.xy_lb_size)
    img.axs[0][i].legend(fontsize=img.legend_size, loc="best", frameon=False)

# img.subimg_adjust(h=0, w=0)
img.save_img(pic_path)
print(pic_path)
img.close_img()


