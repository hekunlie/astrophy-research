import numpy
from astropy.io import fits
import h5py
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
import tool_box
from plot_tool import Image_Plot
import matplotlib.pyplot as plt
from numpy import fft

size = 60
cen = int(size / 2)
my, mx = numpy.mgrid[0:size, 0:size]
ky, kx = my - cen, mx - cen

D_f = (kx ** 2 - ky ** 2 + 2j * ky * kx) / (ky ** 2 + kx ** 2)*numpy.pi
D_f_con = (kx ** 2 - ky ** 2 - 2j * ky * kx) / (ky ** 2 + kx ** 2)*numpy.pi
D = D_f*D_f_con
plt.imshow(D.real)
plt.colorbar()
plt.show()
plt.close()
# def shear2kappa(gamma1, gamma2):
#
#     size = gamma1.shape[0]
#     cen = int(size / 2)
#     my, mx = numpy.mgrid[0:size, 0:size]
#     ky, kx = my - cen, mx - cen
#
#     # D_f = (kx ** 2 - ky ** 2 + 2j * ky * kx) / (ky ** 2 + kx ** 2)*numpy.pi
#     D_f_con = (kx**2 - ky**2 - 2j*ky*kx)/(ky**2 + kx**2)*numpy.pi
#     D_f_con[cen,cen] = 0-0j
#
#     gamma_ = gamma1 + 1j*gamma2
#     gamma_f = fft.fftshift(fft.fft2(gamma_))
#
#     kappa_f_recon = gamma_f * D_f_con / numpy.pi
#     kappa_recon = fft.ifft2(fft.ifftshift(kappa_f_recon))
#     return kappa_recon
#
# def kappa2shear(kappa):
#
#     size = kappa.shape[0]
#     cen = int(size / 2)
#     my, mx = numpy.mgrid[0:size, 0:size]
#     ky, kx = my - cen, mx - cen
#     kappa_f = fft.fftshift(fft.fft2(kappa))
#     D_f = (kx ** 2 - ky ** 2 + 2j * ky * kx) / (ky ** 2 + kx ** 2)*numpy.pi
#     D_f[cen, cen] = 0 + 0j
#     gamma_f = D_f * kappa_f / numpy.pi
#     gamma = fft.ifft2(fft.ifftshift(gamma_f))
#     return gamma


size = 40
cen = int(size/2)
sig = 5
my, mx = numpy.mgrid[0:size,0:size]
ky, kx = my - cen, mx-cen

kappa = 10/numpy.pi/2/sig**2*numpy.exp(-((my-cen)**2+(mx-cen)**2)/2/sig**2)

kappa_f = fft.fftshift(fft.fft2(kappa))
k_pow = numpy.abs(kappa_f)**2

D_f = (kx ** 2 - ky ** 2 + 2j * ky * kx) / (ky ** 2 + kx ** 2)

print(kappa.sum())
print(k_pow.max())

gamma = tool_box.kappa2shear(kappa)

g1, g2 = gamma.real, gamma.imag

kappa_recon = tool_box.shear2kappa(g1, g2).real

img = Image_Plot()
img.create_subfig(1,3)
fig1 = img.axs[0][0].imshow(kappa)
plt.colorbar(fig1, ax=img.axs[0][0])
fig2 = img.axs[0][1].imshow(k_pow)
plt.colorbar(fig2, ax=img.axs[0][1])
fig3 = img.axs[0][2].imshow(numpy.abs(D_f))
plt.colorbar(fig3, ax=img.axs[0][2])
img.show_img()


g = numpy.abs(gamma)
npw = numpy.where(g == 0)

cos_theta = numpy.sqrt((1+g1/g)/2)
sin_theta = g2/2/g/cos_theta
idx = cos_theta == 0
sin_theta[idx] = 1

img = Image_Plot()
img.create_subfig(1,3)
fig1 = img.axs[0][0].imshow(g)
plt.colorbar(fig1, ax=img.axs[0][0])
fig2 = img.axs[0][1].imshow(g1)
plt.colorbar(fig2, ax=img.axs[0][1])
fig3 = img.axs[0][2].imshow(g2)
plt.colorbar(fig3, ax=img.axs[0][2])
img.show_img()

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

shear_bench = 0.02
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

img = Image_Plot()
img.create_subfig(2,2)

fig1 = img.axs[0][0].imshow(kappa)
plt.colorbar(fig1, ax=img.axs[0][0])

fig2 = img.axs[0][1].imshow(kappa_recon)
plt.colorbar(fig2, ax=img.axs[0][1])

fig3 = img.axs[1][0].imshow(kappa_recon-kappa)
plt.colorbar(fig3, ax=img.axs[1][0])

diff = (kappa - kappa_recon)/kappa
fig4 = img.axs[1][1].imshow(diff)
plt.colorbar(fig4, ax=img.axs[1][1])
img.show_img()