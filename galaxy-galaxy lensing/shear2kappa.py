import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
path.append('%s/work/mylib/'%my_home)
import numpy
import tool_box
from plot_tool import Image_Plot
import matplotlib.pyplot as plt
from numpy import fft
from scipy import signal


def gauss_kappa(ny, nx, ceny, cenx, sig, ampl=1):
    cy, cx = int(ny / 2), int(nx / 2)
    ky, kx = numpy.mgrid[-cy:cy, -cx:cx]
    return ampl/numpy.pi/2/sig**2*numpy.exp(-((ky-ceny)**2+(kx-cenx)**2)/2/sig**2)

def add_hole(ny, nx, ceny, cenx, radius):
    my, mx = numpy.mgrid[0:ny, 0:nx]
    mask_final = numpy.zeros_like(my)
    hole_num = len(cenx)
    for i in range(hole_num):
        mask = (my-ceny[i])**2 + (mx-cenx[i])**2
        idx = mask <= radius[i]**2
        mask[idx] = 0
        idx = mask > 0
        mask[idx] = 1
        mask_final += mask
    idx = mask_final < hole_num - 1.e-5
    mask_final[idx] = 0
    return mask_final


def image_fft(image):
    return fft.fftshift(fft.fft2(image))

def image_ifft(image):
    return fft.ifft2(fft.ifftshift(image))

def shear2kappa(gamma1, gamma2, cen_x=0, cen_y=0):
    size = gamma1.shape[0]
    cen = int(size / 2)
    my, mx = numpy.mgrid[0:size, 0:size]
    ky, kx = my - cen, mx - cen

    # D_f = (kx ** 2 - ky ** 2 + 2j * ky * kx) / (ky ** 2 + kx ** 2)*numpy.pi
    D_f_con = (kx**2 - ky**2 - 2j*ky*kx)/(ky**2 + kx**2)*numpy.pi

    D_f_con[cen,cen] = cen_x-cen_y*1j

    gamma_ = gamma1 + 1j*gamma2
    gamma_f = image_fft(gamma_)

    kappa_f_recon = gamma_f * D_f_con / numpy.pi
    kappa_recon = image_ifft(kappa_f_recon)
    return kappa_recon

def kappa2shear(kappa, cen_x=0, cen_y=0):
    size = kappa.shape[0]
    cen = int(size / 2)
    my, mx = numpy.mgrid[0:size, 0:size]
    ky, kx = my - cen, mx - cen

    kappa_f = image_fft(kappa)
    D_f = (kx ** 2 - ky ** 2 + 2j * ky * kx) / (ky ** 2 + kx ** 2)*numpy.pi
    D_f[cen, cen] = cen_x-cen_y*1j

    gamma_f = D_f * kappa_f / numpy.pi
    gamma = image_ifft(gamma_f)
    return gamma

def gamma_to_kappa(g1, g2):
    sp = g1.shape
    nx, ny = sp[1]*2+1, sp[0]*2+1
    cenx, ceny = int(nx / 2), int(ny / 2)

    my, mx = numpy.mgrid[-ceny:ny - ceny, -cenx:nx - cenx]
    print(nx, cenx, my.min(), my.max())
    R2 = my ** 2 + mx ** 2
    D1 = (my ** 2 - mx ** 2) / R2 ** 2
    D2 = -2 * my * mx / R2 ** 2

    D1[ceny, cenx] = 0
    D2[ceny, cenx] = 0

    kappa_1 = signal.convolve(g1, D1, mode="same")
    kappa_2 = signal.convolve(g2, D2, mode="same")

    # kappa_ini = numpy.zeros_like(g1) + 0.001

    # while True:
    #     max_kappa_ini = kappa_ini.max()
    #
    #     g1_prime = (1 - kappa_ini)*g1
    #     g2_prime = (1 - kappa_ini)*g2
    #
    #     kappa_1 = signal.convolve(g1_prime, D1, mode="same")
    #     kappa_2 = signal.convolve(g2_prime, D2, mode="same")
    #
    #     kappa_resc = kappa_1 + kappa_2
    #
    #
    #     kappa_ini = kappa_resc
    #
    #     diff_kappa = kappa_resc.max() - max_kappa

    return (kappa_1 + kappa_2)/numpy.pi

size = 50
cen = int(size/2)
sig = 5


mask = add_hole(size, size, [10, 15, 40], [10, 16, 24], [3, 4, 7])
plt.imshow(mask)
plt.show()
exit()
kappa_1 = gauss_kappa(size, size, 14, -8,sig)
kappa_2 = gauss_kappa(size, size, -10, 2,sig)
kappa = kappa_1 + kappa_2
gamma = kappa2shear(kappa)
g1, g2 = gamma.real, gamma.imag
kappa_recon = shear2kappa(g1,g2).real
kappa_recon_c = gamma_to_kappa(g1, g2)


g = numpy.abs(gamma)
npw = numpy.where(g == 0)

cos_theta = numpy.sqrt((1+g1/g)/2)
sin_theta = g2/2/g/cos_theta
idx = cos_theta == 0
sin_theta[idx] = 1
plot_data_1 = [[kappa, g],[g1,g2]]

img = Image_Plot()
img.create_subfig(2,2)

for i in range(2):
    for j in range(2):
        fig = img.axs[i][j].imshow(plot_data_1[i][j])
        plt.colorbar(fig, ax=img.axs[i][j])
img.show_img()



# plot shear map
img = Image_Plot(fig_x=10, fig_y=10)
img.create_subfig(1, 1)

max_g = g.max()

ra_bin = numpy.linspace(-cen - 1, cen, size + 1)
dec_bin = numpy.linspace(-cen - 1, cen, size + 1)

nx, ny = size, size

ra_min, ra_max = ra_bin.min(), ra_bin.max()
dec_min, dec_max = dec_bin.min(), dec_bin.max()

ra_sq_len = ra_bin[2] - ra_bin[1]
dec_sq_len = ra_sq_len
max_len = ra_sq_len * 0.9

shear_bench = 0.002
scale_len = shear_bench / max_g * max_len
print(max_g)
x1, x2 = ra_min + ra_sq_len * 6, ra_min + ra_sq_len * 6 + scale_len
y1, y2 = dec_max + dec_sq_len, dec_max + dec_sq_len

img.axs[0][0].plot([x1, x2], [y1, y2], c="black")
img.axs[0][0].text(0.04, 0.95, "$\gamma=%.3f$" % shear_bench, color='black', ha='left',
                   va='center', transform=img.axs[0][0].transAxes, fontsize=img.legend_size - 5)

for i in range(ny + 1):
    img.axs[0][0].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--", alpha=0.6, linewidth=0.3)
for j in range(nx + 1):
    img.axs[0][0].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black", linestyle="--", alpha=0.6, linewidth=0.3)

dg_scale = g / max_g * max_len / 2

for i in range(ny):
    for j in range(nx):
        if g2[i, j] < 0:
            dx = -dg_scale[i, j] * cos_theta[i, j]
            dy = -dg_scale[i, j] * sin_theta[i, j]
        else:
            dx = dg_scale[i, j] * cos_theta[i, j]
            dy = dg_scale[i, j] * sin_theta[i, j]

        x = (ra_bin[j] + ra_bin[j + 1]) / 2
        y = (dec_bin[i] + dec_bin[i + 1]) / 2
        img.axs[0][0].plot([x + dx, x - dx], [y + dy, y - dy], c="C0")
img.show_img()


plt_data = [[kappa_recon_c, kappa, kappa_recon],
            [kappa - kappa_recon_c, kappa - kappa, kappa_recon-kappa]]
print(plt_data[1][2].max(), plt_data[1][2].min())
print(plt_data[1][2])
img = Image_Plot(fig_x=8, fig_y=8)
img.create_subfig(2, 3)
for i in range(2):
    for j in range(3):
        subfig = img.axs[i][j].imshow(plt_data[i][j])
        plt.colorbar(subfig, ax=img.axs[i][j])
img.show_img()
