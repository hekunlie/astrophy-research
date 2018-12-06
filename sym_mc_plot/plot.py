import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import time

i = 0.000000001
while True:
    i += 0.00000001
    print(i)
    time.sleep(10)





fig_x = 8
fig_y = fig_x*4/6
figs = (fig_x*2, fig_x*4/3)
fig = plt.figure(figsize=figs)

fonts = 20
xy_lb_size = 16
lenged_size = fonts - 6
axis_linewidth = 1.2


total_path = "F:\works\selection_bias\data\galsim\dimmer/"
total_path = total_path.replace("\\","/")
filter_names = ["sex2_", "sex3_", "sex4_"]
select = ["flux2", "sex_snr", "mag_auto"]
lbs = ["SNR$_F$", "SNR$_S$", "MAG"]
gauss_filter = ["gauss_2.0", "gauss_3.0", "gauss_4.0"]
xfmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(xfmt)
ch = numpy.arange(7)


# dimmer c/m
# xy_lim = [[-0.00005, -0.00005, -0.00005], [0.00005, 0.00005, 0.00005]]
# text_y = [-0.004,-0.009,-0.004]
# y_ticks = [numpy.arange(-0.00004, 0.00005, 0.00002),
#            numpy.arange(-0.00004, 0.00005, 0.00002),
#            numpy.arange(-0.00004, 0.00005, 0.00002)]
xy_lim = [[-0.0045, -0.012, -0.0045], [0.0025, 0.002, 0.0035]]
text_y = [-0.004,-0.009,-0.004]
y_ticks = [numpy.arange(-0.004, 0.0021, 0.002),
           numpy.arange(-0.009, 0.002, 0.003),
           numpy.arange(-0.004, 0.0021, 0.002)]

# dimmerm3 c / m
# xy_lim = [[-0.00006, -0.00006, -0.00006], [0.00005, 0.00005, 0.00005]]
# text_y = [-0.004,-0.009,-0.004]
# y_ticks = [numpy.arange(-0.00004, 0.00005, 0.00002),
#            numpy.arange(-0.00004, 0.00005, 0.00002),
#            numpy.arange(-0.00004, 0.00005, 0.00002)]
# xy_lim = [[-0.0045, -0.013, -0.0055], [0.005, 0.0024, 0.006]]
# text_y = [-0.004,-0.009,-0.004]
# y_ticks = [numpy.arange(-0.004, 0.0041, 0.002),
#            numpy.arange(-0.01, 0.002, 0.0025),
#            numpy.arange(-0.005, 0.0051, 0.0025)]
# pts 64 m / c
# xy_lim = [[-0.023, -0.045, -0.035], [0.012, 0.01, 0.023]]
# text_y = [-0.004,-0.009,-0.004]
# y_ticks = [numpy.arange(-0.02, 0.011, 0.01),
#            numpy.arange(-0.04, 0.01, 0.01),
#            numpy.arange(-0.03, 0.016, 0.015)]
# xy_lim = [[-0.00022, -0.00022, -0.00022], [0.00015, 0.00015, 0.00015]]
# text_y = [-0.004,-0.009, -0.004]
# y_ticks = [numpy.arange(-0.0002, 0.00015, 0.0001),
#            numpy.arange(-0.0002, 0.00015, 0.0001),
#            numpy.arange(-0.0002, 0.00015, 0.0001)]

# pts 64_bright m / c
# xy_lim = [[-0.018, -0.04, -0.018], [0.012, 0.005, 0.012]]
# text_y = [-0.004,-0.009,-0.004]
# y_ticks = [numpy.arange(-0.016, 0.011, 0.008),
#            numpy.arange(-0.03, 0.01, 0.01),
#            numpy.arange(-0.016, 0.011, 0.008)]

# xy_lim = [[-0.00014, -0.00014, -0.00014], [0.00012, 0.00012, 0.00012]]
# text_y = [-0.004,-0.009, -0.004]
# y_ticks = [numpy.arange(-0.0001, 0.00015, 0.00005),
#            numpy.arange(-0.0001, 0.00015, 0.00005),
#            numpy.arange(-0.0001, 0.00015, 0.00005)]

sig = "2"
mc_plot = 'm'
lim = 1
if mc_plot == "m":
    plot_tag = 0
    stand = 1
else:
    plot_tag = 2
    stand = 0

pic_name = total_path+"%s_%ssigma.pdf"%(mc_plot,sig)
markers = ["s", "p", '>']
colors = ["blue", "red", 'green']
for row in range(3):
    if row == 0:
        data_path = total_path + "sex2_%s/flux2/total.npz" % sig
        data = numpy.load(data_path)
        mcs = [data["arr_0"], data["arr_1"]]
        for col in range(2):
            ax = fig.add_subplot(3, 2, 2*row+col+1)
            lb = "SNR$_F$: $%s_%d$"%(mc_plot, col+1)
            ax.errorbar(ch*10.,mcs[col][:,ch][plot_tag] - stand, mcs[col][:,ch][plot_tag+1], color=colors[0],capsize=4,
                        marker='s',ms=5, label=lb)

            x = ax.set_xlim()
            if lim:
                y = ax.set_ylim(xy_lim[0][row], xy_lim[1][row])
                ax.set_yticks(y_ticks[row])
            if col == 0:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=1,fontsize=lenged_size, loc='upper left',edgecolor="black")
            else:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=1,fontsize=lenged_size, loc='upper left',edgecolor="black")
            ax.plot([x[0], x[1]],[0,0], linestyle="-.", c='grey')
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)

            if col == 1:
                ax.set_yticklabels([])
            for axis in ["bottom", "left", "top", "right"]:
                ax.spines[axis].set_linewidth(axis_linewidth)
            ax.xaxis.set_tick_params(which="both", direction="in", length=4, width=axis_linewidth)

    elif row == 1:
        mcs = [[],[]]
        for i in range(3):
            data_path = total_path + filter_names[i] + "%s/%s/total.npz" % (sig,select[row])
            data = numpy.load(data_path)
            mcs[0].append(data['arr_0'])
            mcs[1].append(data['arr_1'])
        for col in range(2):
            ax = fig.add_subplot(3, 2, 2*row+col+1)
            for ii in range(3):
                lb = "SNR$_S$: $%s_%d$, %s" % (mc_plot,col+1, gauss_filter[ii])
                ax.errorbar(ch*10., mcs[col][ii][:,ch][plot_tag] - stand, mcs[col][ii][:,ch][plot_tag+1],color=colors[ii],
                            capsize=4, marker=markers[ii],ms=5,label=lb)
            if col == 0:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=1,fontsize=lenged_size, loc='upper right',edgecolor="black",facecolor="white")#,bbox_to_anchor=(0.98,0.7))
            else:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=1,fontsize=lenged_size, loc='upper right',edgecolor="black",facecolor="white")#,bbox_to_anchor=(0.98,0.7))
            ax.set_xlim(x[0],x[1])
            if lim:
                y = ax.set_ylim(xy_lim[0][row], xy_lim[1][row])
                ax.set_yticks(y_ticks[row])
            ax.plot([x[0], x[1]],[0,0], linestyle="-.", c='grey')
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)

            if col == 1:
                ax.set_yticklabels([])
            for axis in ["bottom", "left", "top", "right"]:
                ax.spines[axis].set_linewidth(axis_linewidth)
            ax.xaxis.set_tick_params(which="both", direction="in", length=4, width=axis_linewidth)
    else:
        mcs = [[],[]]
        for i in range(3):
            data_path = total_path + filter_names[i] + "%s/%s/total.npz" % (sig, select[row])
            print(os.path.exists(data_path))
            data = numpy.load(data_path)
            mcs[0].append(data['arr_0'])
            mcs[1].append(data['arr_1'])
        for col in range(2):
            ax = fig.add_subplot(3, 2, 2*row+col+1)
            for ii in range(3):
                lb = "MAG: $%s_%d$, %s" % (mc_plot,col+1, gauss_filter[ii])
                ax.errorbar(ch * 10., mcs[col][ii][:, ch][plot_tag] - stand, mcs[col][ii][:, ch][plot_tag + 1],
                            color=colors[ii],
                            capsize=4, marker=markers[ii], ms=5, label=lb)
            if col == 0:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc='upper left',edgecolor="black",facecolor="white")
            else:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc='upper left',edgecolor="black",facecolor="white")
            ax.set_xlim(x[0], x[1])
            if lim:
                y = ax.set_ylim(xy_lim[0][row], xy_lim[1][row])
                ax.set_yticks(y_ticks[row])
            ax.plot([x[0], x[1]], [0, 0], linestyle="-.", c='grey')
            ax.xaxis.set_major_formatter(xticks)
            ax.set_xlabel("Cutoff percentage", fontsize=fonts)
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            if col == 1:
                ax.set_yticklabels([])
            for axis in ["bottom", "left", "top", "right"]:
                ax.spines[axis].set_linewidth(axis_linewidth)
            ax.xaxis.set_tick_params(which="both", direction="in", length=4, width=axis_linewidth)

plt.subplots_adjust(wspace=0,hspace=0)
plt.savefig(pic_name,bbox_inches='tight')
plt.show()