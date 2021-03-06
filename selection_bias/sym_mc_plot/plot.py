import numpy
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
import os
import time



fig_x = 16
fig_y = fig_x*4/6
figs = (fig_x*2, fig_x*4/3)
fig = plt.figure(figsize=figs)

fonts = 20
xy_lb_size = 16
lenged_size = fonts - 10
axis_linewidth = 1.2
cap_size = 5
line_w = 2
total_path = "E:\cuts\sym/"
total_path = total_path.replace("\\","/")
filter_names = ["sex2_","sex3_",  "sex4_"]
select = ["flux2", "sex_snr", "mag_auto","snr_auto"]
lbs = ["SNR$_F$", "SNR$_S$", "MAG","SNR$_A$"]
ex_lbs = ["P$_{k0}$", "P$_{k0,fit}$", "MAX(P$_{k0}$, P$_{k0,fit}$","MAX(SNR$_F$, SNR$_F$-fit", "MAG$_{true}$"]
gauss_filter = ["gauss_2.0", "gauss_3.0", "gauss_4.0"]
xfmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(xfmt)
ch = numpy.arange(10)



xy_lim = [[-0.008, -0.015, -0.01], [0.0025, 0.002, 0.0035]]
text_y = [-0.004,-0.009,-0.004]
y_ticks = [numpy.arange(-0.004, 0.0021, 0.002),
           numpy.arange(-0.009, 0.002, 0.003),
           numpy.arange(-0.004, 0.0021, 0.002)]

sig = "1.5"
mc_plot = 'm'
lim = 0
if mc_plot == "m":
    plot_tag = 0
    stand = 1
else:
    plot_tag = 2
    stand = 0

pic_name = total_path+"%s_%ssigma.png"%(mc_plot,sig)
markers = ["s", "p", '>']
colors = ["blue", "red", 'brown']
leg_pos = ['best','best','best','best','best','best','best','best']
fig_rows = 4
for row in range(fig_rows):
    if row == 0:
        data_path = total_path + "sex2_%s/flux2/total.npz" % sig
        data = numpy.load(data_path)
        mcs = [data["arr_0"], data["arr_1"]]
        try:
            data_path = total_path + "sex2_%s/flux_alt/total.npz" % sig

            data = numpy.load(data_path)
            mcs_alt = [data["arr_0"], data["arr_1"]]
            flux_exs_mcs = []
            for i2 in range(5):
                data_path = total_path + "sex2_%s/flux2_ex%d/total.npz" %(sig, i2+1)
                print(os.path.exists(data_path))
                data = numpy.load(data_path)
                flux_exs_mcs.append([data["arr_0"], data["arr_1"]])
        except:
            print("No flux_alt")

        for col in range(2):
            ax = fig.add_subplot(fig_rows, 2, 2*row+col+1)
            lb = "SNR$_F$: $%s_%d$"%(mc_plot, col+1)
            pts_show = mcs[col][:,ch][plot_tag] - stand
            pts_err = mcs[col][:,ch][plot_tag+1]
            ax.errorbar(ch*10.,pts_show, pts_err, color=colors[0],linewidth=line_w,capsize=cap_size, marker='s',ms=5, label=lb)
            try:
                lb = "SNR$_F$-fit: $%s_%d$"%(mc_plot, col+1)
                pts_show = mcs_alt[col][:,ch][plot_tag] - stand
                pts_err = mcs_alt[col][:,ch][plot_tag+1]
                ax.errorbar(ch*10.,pts_show, pts_err, color=colors[2],linewidth=line_w,capsize=cap_size, marker='s',ms=5, label=lb)
                for i2 in range(5):
                    pts_show = flux_exs_mcs[i2][col][:, ch][plot_tag] - stand
                    pts_err = flux_exs_mcs[i2][col][:, ch][plot_tag + 1]
                    ax.errorbar(ch * 10., pts_show, pts_err, color="C%d"%i2, linewidth=line_w,capsize=cap_size, marker='s', ms=5, label=ex_lbs[i2])
            except:
                print("No flux_alt")

            x = ax.set_xlim()
            if lim:
                y = ax.set_ylim(xy_lim[0][row], xy_lim[1][row])
                # ax.set_yticks(y_ticks[row])
            if col == 0:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=3,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black")
            else:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=3,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black")
            ax.plot([x[0], x[1]],[0,0], linestyle="-.", c='grey')
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            ax.xaxis.set_major_formatter(xticks)
            #
            # if col == 1:
            #     ax.set_yticklabels([])
            ax.set_xticklabels([])
            for axis in ["bottom", "left", "top", "right"]:
                ax.spines[axis].set_linewidth(axis_linewidth)
            ax.xaxis.set_tick_params(which="both", direction="in", length=cap_size, width=axis_linewidth)

    elif row == 1:
        mcs = [[],[]]
        for i in range(len(filter_names)):
            data_path = total_path + filter_names[i] + "%s/%s/total.npz" % (sig,select[row])
            data = numpy.load(data_path)
            mcs[0].append(data['arr_0'])
            mcs[1].append(data['arr_1'])
        for col in range(2):
            ax = fig.add_subplot(fig_rows, 2, 2*row+col+1)
            for ii in range(len(filter_names)):
                lb = "SNR$_S$: $%s_%d$, %s" % (mc_plot,col+1, gauss_filter[ii])
                pts_show = mcs[col][ii][:,ch][plot_tag] - stand
                pts_err = mcs[col][ii][:,ch][plot_tag+1]
                ax.errorbar(ch*10., pts_show, pts_err,color=colors[ii],linewidth=line_w,capsize=cap_size, marker=markers[ii],ms=5,label=lb)
            if col == 0:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black",facecolor="white")#,bbox_to_anchor=(0.98,0.7))
            else:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black",facecolor="white")#,bbox_to_anchor=(0.98,0.7))
            ax.set_xlim(x[0],x[1])
            if lim:
                y = ax.set_ylim(xy_lim[0][row], xy_lim[1][row])
                # ax.set_yticks(y_ticks[row])
            ax.plot([x[0], x[1]],[0,0], linestyle="-.", c='grey')
            ax.xaxis.set_major_formatter(xticks)
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            ax.set_xticklabels([])
            # if col == 1:
            #     ax.set_yticklabels([])
            for axis in ["bottom", "left", "top", "right"]:
                ax.spines[axis].set_linewidth(axis_linewidth)
            ax.xaxis.set_tick_params(which="both", direction="in", length=4, width=axis_linewidth)
    elif row == 2:
        mcs = [[],[]]
        for i in range(len(filter_names)):
            data_path = total_path + filter_names[i] + "%s/%s/total.npz" % (sig, select[row])
            data = numpy.load(data_path)
            mcs[0].append(data['arr_0'])
            mcs[1].append(data['arr_1'])
        for col in range(2):
            ax = fig.add_subplot(fig_rows, 2, 2*row+col+1)
            for ii in range(len(filter_names)):
                lb = "MAG: $%s_%d$, %s" % (mc_plot,col+1, gauss_filter[ii])
                pts_show = mcs[col][ii][:, ch][plot_tag] - stand
                pts_err = mcs[col][ii][:, ch][plot_tag + 1]
                ax.errorbar(ch * 10., pts_show,pts_err, color=colors[ii], linewidth=line_w,capsize=cap_size, marker=markers[ii], ms=5, label=lb)
            if col == 0:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black",facecolor="white")
            else:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black",facecolor="white")
            ax.set_xlim(x[0], x[1])
            if lim:
                y = ax.set_ylim(xy_lim[0][row], xy_lim[1][row])
                # ax.set_yticks(y_ticks[row])
            ax.plot([x[0], x[1]], [0, 0], linestyle="-.", c='grey')
            ax.xaxis.set_major_formatter(xticks)
            # ax.set_xlabel("Cutoff percentage", fontsize=fonts)
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            # if col == 1:
            #     ax.set_yticklabels([])
            for axis in ["bottom", "left", "top", "right"]:
                ax.spines[axis].set_linewidth(axis_linewidth)
            ax.xaxis.set_tick_params(which="both", direction="in", length=cap_size, width=axis_linewidth)
            ax.set_xticklabels([])
    elif row == 3:
        mcs = [[],[]]
        for i in range(len(filter_names)):
            data_path = total_path + filter_names[i] + "%s/%s/total.npz" % (sig, select[row])
            data = numpy.load(data_path)
            mcs[0].append(data['arr_0'])
            mcs[1].append(data['arr_1'])
        for col in range(2):
            ax = fig.add_subplot(fig_rows, 2, 2*row+col+1)
            for ii in range(len(filter_names)):
                lb = "SNR$_A$: $%s_%d$, %s" % (mc_plot,col+1, gauss_filter[ii])
                pts_show = mcs[col][ii][:, ch][plot_tag] - stand
                pts_err = mcs[col][ii][:, ch][plot_tag + 1]
                ax.errorbar(ch * 10., pts_show,pts_err, color=colors[ii], linewidth=line_w,capsize=cap_size, marker=markers[ii], ms=5, label=lb)
            if col == 0:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black",facecolor="white")
            else:
                # ax.text(0, text_y[row],lbs[row]+": $%s_%d$"%(mc_plot, col+1),fontsize=fonts-4)
                ax.legend(ncol=2,fontsize=lenged_size, loc=leg_pos[row*2+col],edgecolor="black",facecolor="white")
            ax.set_xlim(x[0], x[1])
            if lim:
                y = ax.set_ylim(xy_lim[0][row], xy_lim[1][row])
                # ax.set_yticks(y_ticks[row])
            ax.plot([x[0], x[1]], [0, 0], linestyle="-.", c='grey')
            ax.xaxis.set_major_formatter(xticks)
            ax.set_xlabel("Cutoff percentage", fontsize=fonts)
            ax.tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
            # if col == 1:
            #     ax.set_yticklabels([])
            for axis in ["bottom", "left", "top", "right"]:
                ax.spines[axis].set_linewidth(axis_linewidth)
            ax.xaxis.set_tick_params(which="both", direction="in", length=4, width=axis_linewidth)

plt.subplots_adjust(wspace=0.1,hspace=0.1)
plt.savefig(pic_name,bbox_inches='tight')
plt.show()