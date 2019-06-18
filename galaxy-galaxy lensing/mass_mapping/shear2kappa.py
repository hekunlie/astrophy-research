import platform
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append("E:/Github/astrophy-research/mylib/")
path.append('%s/work/mylib/'%my_home)
import numpy
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


img_path = "./pic/"
size = int(argv[1])
sig_level = float(argv[2])
cen = int(size/2)
sig = 5
inverse = range(size-1, -1, -1)

# mask = add_hole(size, size, [10, 15, 40], [10, 16, 24], [3, 4, 7])
mask = add_hole(size, size, [cen-14 - 8], [cen-12-8], [5])

# plt.imshow(mask)
# plt.show()
# exit()

# create kappa
kappa_1 = gauss_kappa(size, size, -14, -12, sig, 10)
kappa_2 = gauss_kappa(size, size, 7, 8, sig, 10)
kappa = (kappa_1 + kappa_2)
kappa = gauss_kappa(size, size, 0,0, sig, 10)

# kappa to gamma = g1 + i*g2
gamma = kappa2shear(kappa)
g1, g2 = gamma.real, gamma.imag
# add noise to each component
g1_noise = numpy.random.normal(0, numpy.abs(g1)/sig_level)
g2_noise = numpy.random.normal(0, numpy.abs(g2)/sig_level)
g1_n = g1 + g1_noise
g2_n = g2 + g2_noise
gamma_n = g1_n + 1j*g2_n

# inverse (KS) method, gamma to kappa
kappa_recon = shear2kappa(g1, g2).real
kappa_recon_n = shear2kappa(g1_n, g2_n).real

# reconstruct kappa field with gamma
kappa_recon_c = gamma_to_kappa(g1, g2)
kappa_recon_c_n = gamma_to_kappa(g1_n, g2_n)


g = numpy.abs(gamma)
npw = numpy.where(g == 0)
cos_theta = numpy.sqrt((1+g1/g)/2)
sin_theta = g2/2/g/cos_theta
idx = cos_theta == 0
sin_theta[idx] = 1

g_n = numpy.abs(gamma_n)
npw_n = numpy.where(g_n == 0)
cos_theta_n = numpy.sqrt((1+g1_n/g_n)/2)
sin_theta_n = g2_n/2/g_n/cos_theta_n
idx_n = cos_theta_n == 0
sin_theta_n[idx] = 1

plot_data_1 = [[kappa, g, g1, g2],
               [kappa_recon, kappa_recon - kappa, kappa_recon_c,kappa_recon_c - kappa]]
plot_data_2 = [[kappa, g_n, g1_n, g2_n],
               [kappa_recon_n, kappa_recon_n - kappa, kappa_recon_c_n, kappa_recon_c_n - kappa]]
titles_1 = [["$\kappa$", "$\gamma$", "$\gamma_1$", "$\gamma_2$"],
            ["$\kappa$ Fourerier", "residual %.6f"%(kappa_recon - kappa).max(), "$\kappa$ real space", "residual"]]
titles_2 = [["$\kappa$", "$\gamma$", "$\gamma_1$", "$\gamma_2$"],
            ["$\kappa$ Fourerier", "residual", "$\kappa$ real space", "residual"]]

img = Image_Plot()
img.subplots(2,4)

for i in range(2):
    for j in range(4):
        fig = img.axs[i][j].imshow(plot_data_1[i][j][inverse])
        img.axs[i][j].set_title(titles_1[i][j], fontsize=img.xy_lb_size)
        plt.colorbar(fig, ax=img.axs[i][j])
img.save_img(img_path + "noise_free.png")
if platform.system() != 'Linux':
    img.show_img()
img.close_img()

img = Image_Plot()
img.subplots(2,4)

for i in range(2):
    for j in range(4):
        fig = img.axs[i][j].imshow(plot_data_2[i][j][inverse])
        img.axs[i][j].set_title(titles_2[i][j], fontsize=img.xy_lb_size)
        plt.colorbar(fig, ax=img.axs[i][j])
img.save_img(img_path + "noisy.png")
if platform.system() != 'Linux':
    img.show_img()
img.close_img()

# plot shear map
img = Image_Plot(fig_x=15, fig_y=15)
img.subplots(1, 2)
img.axs[0][0].set_title("Input", fontsize=img.xy_lb_size)
img.axs[0][1].set_title("Noisy", fontsize=img.xy_lb_size)

ra_bin = numpy.linspace(-cen - 1, cen, size + 1)
dec_bin = numpy.linspace(-cen - 1, cen, size + 1)

nx, ny = size, size

ra_min, ra_max = ra_bin.min(), ra_bin.max()
dec_min, dec_max = dec_bin.min(), dec_bin.max()

ra_sq_len = ra_bin[2] - ra_bin[1]
dec_sq_len = ra_sq_len
max_len = ra_sq_len * 0.9

max_g = g.max()
max_g_n = g_n.max()

shear_bench = 0.02
scale_len = shear_bench / max_g * max_len
scale_len_n = shear_bench / max_g_n * max_len
print(max_g, max_g_n)

x1, x2 = ra_min + ra_sq_len * 6, ra_min + ra_sq_len * 6 + scale_len
x1_n, x2_n = ra_min + ra_sq_len * 6, ra_min + ra_sq_len * 6 + scale_len
y1, y2 = dec_max + dec_sq_len, dec_max + dec_sq_len

img.axs[0][0].plot([x1, x2], [y1, y2], c="black")
img.axs[0][0].text(0.04, 0.95, "$\gamma=%.3f$" % shear_bench, color='black', ha='left',
                   va='center', transform=img.axs[0][0].transAxes, fontsize=img.legend_size - 5)

for i in range(ny + 1):
    img.axs[0][0].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--", alpha=0.6, linewidth=0.3)
    img.axs[0][1].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--", alpha=0.6, linewidth=0.3)
for j in range(nx + 1):
    img.axs[0][0].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black", linestyle="--", alpha=0.6, linewidth=0.3)
    img.axs[0][1].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black", linestyle="--", alpha=0.6, linewidth=0.3)


dg_scale = g / max_g * max_len / 2
dg_scale_n = g_n / max_g_n *max_len / 2

for i in range(ny):
    for j in range(nx):
        if mask[i,j] > 0:
            if g2[i, j] < 0:
                dx = -dg_scale[i, j] * cos_theta[i, j]
                dy = -dg_scale[i, j] * sin_theta[i, j]

                dx_n = -dg_scale_n[i, j] * cos_theta_n[i, j]
                dy_n = -dg_scale_n[i, j] * sin_theta_n[i, j]
            else:
                dx = dg_scale[i, j] * cos_theta[i, j]
                dy = dg_scale[i, j] * sin_theta[i, j]

                dx_n = dg_scale_n[i, j] * cos_theta_n[i, j]
                dy_n = dg_scale_n[i, j] * sin_theta_n[i, j]

            x = (ra_bin[j] + ra_bin[j + 1]) / 2
            y = (dec_bin[i] + dec_bin[i + 1]) / 2
            img.axs[0][0].plot([x + dx, x - dx], [y + dy, y - dy], c="C0")

            img.axs[0][1].plot([x + dx_n, x - dx_n], [y + dy_n, y - dy_n], c="C0")
img.save_img(img_path + "shear_map.png")
if platform.system() != 'Linux':
    img.show_img()
img.close_img()

