import matplotlib
matplotlib.use("Agg")
from sys import path
path.append('/home/hkli/work/fourier_quad/')
from astropy.io import fits
import numpy
import matplotlib.pyplot as plt
from Fourier_Quad import Fourier_Quad
import tool_box
import lsstetc
from subprocess import Popen
from sys import argv


size, num, title = int(argv[1]), int(argv[2]), argv[3]
seed = 80000

markers = ['o','v','p','h','d','s']
colors = ['red','orange','green','violet','cyan','b']

prop = lsstetc.ETC(band='r', pixel_scale=0.2, stamp_size=size, nvisits=180)
flux = numpy.array([prop.flux(21), prop.flux(24),  prop.flux(24.9)])

fq = Fourier_Quad(size, seed)

fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
ax2 = fig.add_subplot(122)

total_path = '/home/hkli/work/sex/imgs/'

data = numpy.zeros((6, num))
data_0 = numpy.zeros((6, 1))

for i in range(3):
    gal0_path = total_path+'gal0_%d.fits'%i
    gal0 = fits.open(gal0_path)[0].data
    snr_0 = fq.snr_f(gal0)

    fits_name0 = total_path + "gal0_%d.fits" % i
    cat_name0 = total_path + "gal0_%d.cat" % i

    cmd = "sex %s -CATALOG_NAME %s" % (fits_name0, cat_name0)
    sex = Popen(cmd, shell=True)
    sex.wait()

    cata_data0 = numpy.loadtxt(cat_name0)
    cat_row = cata_data0.shape
    sex_snr_0 = cata_data0[0]

    data_0[i*2, 0], data_0[i*2+1, 0] = sex_snr_0, snr_0

    for k in range(num):
        gals_path = total_path + 'gal_%d_%d.fits' %(i, k)
        gal_img = fits.open(gals_path)[0].data
        fits_name = total_path + "gal_%d_%d.fits" %(i, k)
        cat_name = total_path + "gal_%d_%d.cat" %(i, k)
        cmd = "sex %s -CATALOG_NAME %s" % (fits_name, cat_name)
        sex = Popen(cmd, shell=True)
        sex.wait()

        cata_data = numpy.loadtxt(cat_name)
        sex_snr = cata_data[0]
        snr = fq.snr_f(gal_img)
        data[i*2, k], data[i*2+1, k] = sex_snr, snr
    sex_snr_i, snr_i = data[i*2], data[i*2+1]
    sex_delta_i = (sex_snr_i - sex_snr_0)/sex_snr_0
    snr_delta_i = (snr_i - snr_0)/snr_0

    lb1 = "F-SNR=%.4f,S-SNR=%.4f" %(snr_0,sex_snr_0)
    ax1.plot(numpy.linspace(0, numpy.pi, num), sex_delta_i, c=colors[i], ms=6, label=lb1, marker=markers[i],fillstyle='none',linestyle=' ')
    ax2.plot(numpy.linspace(0, numpy.pi, num), snr_delta_i, c=colors[i], ms=6, label=lb1, marker=markers[i],fillstyle='none',linestyle=' ')

y_ticks = numpy.linspace(-10**(-4), 10**(-4), 5)
# ax.set_yticks(y_ticks)
x_ticks = numpy.linspace(0, numpy.pi, 5)
x_ticks_ = [r"0",r"$\frac{\pi}{4}$",r"$\frac{\pi}{2}$",r"$\frac{3\pi}{4}$",r"$\pi$"]
# x_ticks_ = ["0","1","2","3","4"]
ax1.set_xticks(x_ticks)
ax1.set_xticklabels(x_ticks_)
ax2.set_xticks(x_ticks)
ax2.set_xticklabels(x_ticks_)

ax1.tick_params(direction='in', top=True, right=True)
ax2.tick_params(direction='in', top=True, right=True)

ax1.set_xlabel("Theta")
ax2.set_xlabel("Theta")

ax1.yaxis.get_major_formatter().set_powerlimits((1, 2))
ax2.yaxis.get_major_formatter().set_powerlimits((1, 2))


ax1.set_ylabel("Change rate of sex_snr")
ax2.set_ylabel("Change rate of fsnr")
ax1.legend()
ax2.legend()
fig.suptitle(title)

pic_name = total_path + '%s.png'%argv[3]
plt.savefig(pic_name)
plt.show()
plt.close()

data_path = total_path + "total.npz"
numpy.savez(data_path, data_0, data)

