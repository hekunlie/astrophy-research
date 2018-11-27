import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick


fig_x = 7
fig_y = fig_x*4/6
figs = (fig_x*2, fig_x*4/3)
fig = plt.figure(figsize=figs)
fonts = 16
xy_lb_size = 12
total_path = "E:/works/selection_bias/data_for_paper/galsim/dimmer_m3/"
filter_names = ["sex2_", "sex3_", "sex4_"]
select = ["flux2", "sex_snr", "mag_auto"]
lbs = ["F-SNR", "S-SNR", "MAG_AUTO"]
gauss_filter = ["gauss_2.0", "gauss_3.0", "gauss_4.0"]
xfmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(xfmt)
ch = numpy.arange(7)


# for 'm' of /pts/64
# xy_lim = (-0.035, 0.012)
# text_y = -0.03
# for 'c' of /pts/64
# xy_lim = (-0.00024, 0.00022)
# text_y = -0.0002

# for 'm' of /pts/64_bright
# xy_lim = (-0.035, 0.012)
# text_y = -0.02
# for 'c' of /pts/64_bright
# xy_lim = (-0.00024, 0.00015)
# text_y = -0.00015

# for 'm' of /pts/64_dimmer
# xy_lim = (-0.035, 0.012)
# text_y = -0.02
# for 'c' of /pts/64_bright
# xy_lim = (-0.00026, 0.00026)
# text_y = -0.00015


# for 'm' of m1
# xy_lim = (-0.025, 0.012)
# text_y = -0.01
# for 'c' of m3
xy_lim = (-0.0001, 0.00014)
text_y = -0.00007

# for 'c' of m3
# xy_lim = (-0.00022, 0.00022)
# text_y = -0.00013
# for 'c' of m3
# xy_lim = (-0.032, 0.005)
# text_y = -0.022
# for 'c' of m3
# xy_lim = (-0.00022, 0.00022)
# text_y = -0.00013
sig = "2"
mc_plot = 'c'
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
        data = numpy.load(total_path + "flux2_%ssigma/total.npz" % sig)
        mcs = [data["arr_0"], data["arr_1"]]
        for col in range(2):
            ax = fig.add_subplot(321+row*2+col)
            ax.errorbar(ch*10.,mcs[col][:,ch][plot_tag] - stand, mcs[col][:,ch][plot_tag+1], color=colors[0],capsize=4,
                        marker='s',ms=5)
            ax.text(0, text_y, lbs[row] + ": $%s_%d$" % (mc_plot,col + 1), fontsize=fonts - 4)
            x = ax.set_xlim()
            y = ax.set_ylim(xy_lim[0], xy_lim[1])
            ax.plot([x[0], x[1]],[0,0], linestyle="-.", c='black')
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            # ax.set_xticks([])
            if col == 1:
                ax.set_yticklabels([])
    else:
        mcs = [[], []]
        for i in range(3):
            data = numpy.load(total_path + filter_names[i] + "%s/%s/total.npz" % (sig,select[row]))
            mcs[0].append(data['arr_0'])
            mcs[1].append(data['arr_1'])
        for col in range(2):
            ax = fig.add_subplot(321 + row * 2 + col)
            for ii in range(3):
                lb = "%s" % (gauss_filter[ii])
                ax.errorbar(ch*10., mcs[col][ii][:,ch][plot_tag] - stand, mcs[col][ii][:,ch][plot_tag+1],color=colors[ii],
                            capsize=4, marker=markers[ii],ms=5,label=lb)
            ax.text(0, text_y,lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
            if row == 1:
                ax.legend(ncol=3,fontsize=fonts-4, loc='upper left')
            else:
                ax.legend(ncol=3, fontsize=fonts - 4, loc='upper left')
            ax.set_xlim(x[0],x[1])
            y = ax.set_ylim(xy_lim[0], xy_lim[1])
            ax.plot([x[0], x[1]],[0,0], linestyle="-.", c='black')
            if row == 1:
                ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            else:
                ax.xaxis.set_major_formatter(xticks)
                ax.set_xlabel("Cutoff percentage", fontsize=fonts)
                ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            if col == 1:
                ax.set_yticklabels([])


plt.subplots_adjust(wspace=0,hspace=0)
plt.savefig(pic_name)
plt.show()