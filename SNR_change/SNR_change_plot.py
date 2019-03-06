import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from Fourier_Quad import Fourier_Quad
from subprocess import Popen
import tool_box


size, num, title, flux_num = int(argv[1]), int(argv[2]), argv[3], int(argv[4])
seed = 80000

markers = ['o','v','p','h','d','s',"4","*","X","^",">","+"]
colors = ["C%d"%i for i in range(10)]

fmt='%2.f%%'
fig_x = 8
fig_y = fig_x*4/6
figs = (fig_x*4, fig_y)
fonts = 20
xy_lb_size = 18
legend_size = fonts - 4
axis_linewidth = 1.2
plt_line_width = 2
cap_size = 5
fig = plt.figure(figsize=figs)
axs = []
for i in range(4):
    ax = fig.add_subplot(141+i)
    axs.append(ax)
fq = Fourier_Quad(size, seed)
noise_sig = 60
detect_thresh = 1.5
total_path = os.getcwd() + '/imgs/'

data = numpy.zeros((flux_num*4, num))
data_0 = numpy.zeros((flux_num*4, 1))

rim = fq.border(1)
for i in range(flux_num):
    gal0_path = total_path+'gal0_%d.fits'%i
    gal0 = fits.open(gal0_path)[0].data
    snr_0 = numpy.sqrt(fq.pow_spec(gal0)[int(size/2), int(size/2)])/size/noise_sig

    detect_0 = tool_box.stamp_detector(gal0, size, size, 5, 5.5, detect_thresh*noise_sig)
    if detect_0:

        snr_tradi_0 = numpy.sum(detect_0[0]*gal0)/numpy.sqrt(detect_0[0].sum())/noise_sig

        fits_name0 = total_path + "gal0_%d.fits" % i
        cat_name0 = total_path + "gal0_%d.cat" % i
        cmd = "sex %s -CATALOG_NAME %s" % (fits_name0, cat_name0)
        sex = Popen(cmd, shell=True)
        sex.wait()
        try:
            cata_data0 = numpy.loadtxt(cat_name0)
            cat_row = cata_data0.shape
            sex_snr_0 = cata_data0[0]
            flux_auto_0 = cata_data0[1]
            flux_err_0 = cata_data0[2]
            snr_auto_0 = flux_auto_0 / flux_err_0
            sex_mag_0 = cata_data0[3]
        except:
            sex_snr_0, snr_auto_0, sex_mag_0 = 0, 0, 0

        data_0[i*4, 0], data_0[i*4+1, 0], data_0[i*4+2, 0], data_0[i*4+3,0] = sex_snr_0, sex_mag_0, snr_auto_0, snr_0

        for k in range(num):
            gals_path = total_path + 'gal_%d_%d.fits' %(i, k)
            gal_img = fits.open(gals_path)[0].data

            fits_name = total_path + "gal_%d_%d.fits" %(i, k)
            cat_name = total_path + "gal_%d_%d.cat" %(i, k)
            cmd = "sex %s -CATALOG_NAME %s" % (fits_name, cat_name)
            sex = Popen(cmd, shell=True)
            sex.wait()

            try:
                cata_data = numpy.loadtxt(cat_name)
                sex_snr = cata_data[0]
                flux_auto = cata_data[1]
                flux_err = cata_data[2]
                snr_auto = flux_auto / flux_err
                sex_mag = cata_data[3]
            except:
                sex_snr, snr_auto, sex_mag = 0, 0, 0

            snr = numpy.sqrt(fq.pow_spec(gal_img)[int(size / 2), int(size / 2)])/size/noise_sig

            data[i*4, k], data[i*4+1, k], data[i*4+2, k], data[i*4+3,k] = sex_snr, sex_mag, snr_auto, snr

        if sex_mag_0 != 0:

            idx = data[i*4] > 0
            sex_snr_delta = (data[i*4] - sex_snr_0)/sex_snr_0
            sex_mag_delta = (data[i*4+1] - sex_mag_0)/sex_mag_0
            snr_auto_delta = (data[i*4+2] - snr_auto_0) / snr_auto_0
            snr_delta = (data[i*4+3] - snr_0) / snr_0

            deltas = [snr_delta, sex_mag_delta, sex_snr_delta, snr_auto_delta]
            lbs = ["P$_{k0}$", "MAG_AUTO", "SNR$_S$", "SNR$_A$"]
            for select in range(4):
                lb = "%s (%.2f)"%(lbs[select], snr_tradi_0)
                axs[select].plot(numpy.linspace(-0.06, 0.06, num)[idx], deltas[select][idx], c=colors[i], ms=8, label=lb,
                         marker=markers[i],linestyle=' ')#fillstyle='none',

ys = [0,0]
for i in range(4):
    ys = axs[i].set_ylim()
    if ys[1] > ys[1]:
        ys[1] = ys[1]
    if ys[0] < ys[0]:
        ys[0] = ys[0]
dy = ys[1] - ys[0]
x_ticks = numpy.linspace(-0.06, 0.06, 5)
for i in range(4):
    axs[i].set_ylim(-0.042, 0.042)
    axs[i].set_xlim(-0.075,0.075)
    axs[i].set_xticks(x_ticks)
    axs[i].set_xlabel("g1",fontsize=xy_lb_size)
    if i == 0:
        axs[i].set_ylabel("Change rate",fontsize=xy_lb_size)
    else:
        axs[i].set_yticklabels([])
    axs[i].legend(fontsize=legend_size,loc="best")
    axs[i].tick_params(direction='in', labelsize=xy_lb_size, top=True, right=True)
    for axis in ["bottom", "left", "top", "right"]:
        # the line width of the frame
        axs[i].spines[axis].set_linewidth(axis_linewidth)
    axs[i].xaxis.set_tick_params(which="both", direction="in", length=5, width=axis_linewidth)
    axs[i].yaxis.set_tick_params(which="major", direction="in", length=5, width=axis_linewidth)
    axs[i].yaxis.set_tick_params(which="minor", direction="in", length=5, width=axis_linewidth)
plt.subplots_adjust(wspace=0, hspace=0)
pic_name = total_path + '/imgs/%s.pdf'%title
plt.savefig(pic_name,bbox_inches='tight')
plt.show()
plt.close()
print(title)
data_path = total_path + "/imgs/total.npz"
numpy.savez(data_path, data_0, data)

