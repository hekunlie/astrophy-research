import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/' % my_home)
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from Fourier_Quad import Fourier_Quad
from subprocess import Popen
import tool_box

size, num, title, flux_num, step = int(argv[1]), int(argv[2]), argv[3], int(argv[4]), int(argv[5])
seed = 80000

markers = ['o','v','p','h','d','s',"4","*","X","^",">","+"]
colors = ["C%d"%i for i in range(15)]

fq = Fourier_Quad(size, seed)
noise_sig = 60
detect_thresh = 2
fig = plt.figure(figsize=(18,18))
ax1 = fig.add_subplot(221)
ax2 = fig.add_subplot(222)
ax3 = fig.add_subplot(223)
ax4 = fig.add_subplot(224)
total_path = os.getcwd() + '/imgs/'

data = numpy.zeros((flux_num*4, num))
data_0 = numpy.zeros((flux_num*4, 1))

rim = fq.border(1)
for i in range(flux_num):
    gal0_path = total_path+'gal0_%d.fits'%i
    gal0 = fits.open(gal0_path)[0].data
    snr_0 = fq.snr_f(gal0)

    detect_0 = tool_box.stamp_detector(gal0, size, size, 5, 5.5, detect_thresh*noise_sig)
    if detect_0:
        mask, obj_xy = detect_0
        area_0 = len(obj_xy)
        snr_tradition_0 = numpy.sum(gal0*mask)/numpy.sqrt(area_0)/noise_sig
        gal0_pow = fq.pow_spec(gal0)
        noise = numpy.sqrt(numpy.sum(gal0_pow*rim))
        mask_ext = tool_box.edge_extend(mask, size, size, obj_xy, step)

        # plt.imshow(mask_ext)
        # plt.savefig("/home/hkli/work/sex_2/imgs/mask/mask%d_0.png"%i)
        # plt.close()
        idx = mask_ext > 0
        mask_ext[idx] = 1
        snr_0_ext = numpy.abs(numpy.sum(gal0*mask_ext))/noise_sig
    else:
        snr_0_ext = -1
    fits_name0 = total_path + "gal0_%d.fits" % i
    cat_name0 = total_path + "gal0_%d.cat" % i

    cmd = "sex %s -CATALOG_NAME %s" % (fits_name0, cat_name0)
    sex = Popen(cmd, shell=True)
    sex.wait()
    try:
        cata_data0 = numpy.loadtxt(cat_name0)
        cat_row = cata_data0.shape
        sex_snr_0 = cata_data0[0]
        sex_mag_0 = cata_data0[3]
    except:
        sex_snr_0, sex_mag_0 = 0, 0

    data_0[i*4, 0], data_0[i*4+1, 0], data_0[i*4+2, 0], data_0[i*4+3,0] = sex_snr_0, sex_mag_0, snr_0, snr_0_ext

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
            sex_mag = cata_data[3]
        except:
            sex_snr = 0
            sex_mag = 0
        snr = fq.snr_f(gal_img)

        detect = tool_box.stamp_detector(gal_img, size, size, 5, 5.5, detect_thresh*noise_sig)
        if detect:
            mask, obj_xy = detect
            gal_pow = fq.pow_spec(gal_img)
            noise = numpy.sqrt(numpy.sum(gal_pow * rim))
            mask_ext = tool_box.edge_extend(mask, size, size, obj_xy, step)

            # plt.imshow(mask_ext)
            # plt.savefig("/home/hkli/work/sex_2/imgs/mask/mask%d_%d.png"%(i,k))
            # plt.close()
            idx = mask_ext > 0
            mask_ext[idx] = 1
            snr_ext = numpy.abs(numpy.sum(gal_img * mask_ext)) / noise_sig
        else:
            snr_ext = 0

        data[i*4, k], data[i*4+1, k], data[i*4+2, k], data[i*4+3,k] = sex_snr,sex_mag, snr, snr_ext
    sex_snr_i, sex_mag_i, snr_i, snr_ext_i = data[i*4], data[i*4+1], data[i*4+2], data[i*4+3]
    if sex_mag_0 != 0:
        idx = sex_snr_i > 0
        sex_delta_i = (sex_snr_i - sex_snr_0)/sex_snr_0
        snr_delta_i = (snr_i - snr_0)/snr_0
        snr_ext_delta_i = (snr_ext_i - snr_0_ext)/snr_0_ext
        mag_delta_i = (sex_mag_i - sex_mag_0)/sex_mag_0
        idx_s = snr_ext_i > 0
        lb1 = "S-SNR=%.4f" % sex_snr_0
        lb2 = "mag_auto=%.4f" % sex_mag_0
        lb3 = "F-SNR=%.4f" % snr_0
        lb4 = "F-SNR=%.4f" % snr_0_ext
        lb = "SNR: %.3f"%snr_tradition_0
        ax1.plot(numpy.linspace(-0.06, 0.06, num)[idx], sex_delta_i[idx], c=colors[i], ms=6, label=lb,
                 marker=markers[i],fillstyle='none',linestyle=' ')
        ax2.plot(numpy.linspace(-0.06, 0.06, num)[idx], mag_delta_i[idx], c=colors[i], ms=6, label=lb,
                 marker=markers[i],fillstyle='none',linestyle=' ')
        ax3.plot(numpy.linspace(-0.06, 0.06, num)[idx],snr_delta_i[idx], c=colors[i], ms=6, label=lb,
                 marker=markers[i],fillstyle='none',linestyle=' ')
        ax4.plot(numpy.linspace(-0.06, 0.06, num)[idx&idx_s], snr_ext_delta_i[idx&idx_s], c=colors[i], ms=6, label=lb,
                 marker=markers[i], fillstyle='none', linestyle=' ')
y_ticks = numpy.linspace(-10**(-4), 10**(-4), 5)
# ax.set_yticks(y_ticks)
x_ticks = numpy.linspace(-0.06, 0.06, 5)
# x_ticks_ = [r"0",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"]
# x_ticks_ = ["0","1","2","3","4"]
ax1.set_xticks(x_ticks)
# ax1.set_xticklabels(x_ticks_)
ax2.set_xticks(x_ticks)
# ax2.set_xticklabels(x_ticks_)
ax3.set_xticks(x_ticks)
ax4.set_xticks(x_ticks)
# ax4.set_ylim(-0.02,0.02)
ax1.tick_params(direction='in', top=True, right=True)
ax2.tick_params(direction='in', top=True, right=True)
ax3.tick_params(direction='in', top=True, right=True)
ax4.tick_params(direction='in', top=True, right=True)
ax1.set_xlabel("g1")
ax2.set_xlabel("g1")
ax3.set_xlabel("g1")
ax4.set_xlabel("g1")

ax1.yaxis.get_major_formatter().set_powerlimits((1, 2))
ax2.yaxis.get_major_formatter().set_powerlimits((1, 2))
ax3.yaxis.get_major_formatter().set_powerlimits((1, 2))
ax4.yaxis.get_major_formatter().set_powerlimits((1, 2))

ax1.set_ylabel("Change rate of sex_snr")
ax2.set_ylabel("Change rate of mag_auto")
ax3.set_ylabel("Change rate of fsnr")
ax4.set_ylabel("Change rate of fsnr_ext")
ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
fig.suptitle(title)

pic_name = total_path + '/imgs/%s.png'%title
plt.savefig(pic_name,bbox_inches='tight')
plt.show()
plt.close()
print(title)
data_path = total_path + "/imgs/total.npz"
numpy.savez(data_path, data_0, data)

