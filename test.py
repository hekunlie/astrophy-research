# import matplotlib
# matplotlib.use("Agg")
import numpy
import os
# my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
# path.append('%s/work/fourier_quad/'%my_home)
path.append("D:/Github/astrophy-research/mylib/")
import time
from Fourier_Quad import Fourier_Quad
# # import galsim
import matplotlib.pyplot as plt
from astropy.io import fits
import tool_box
import plot_tool
# from mpi4py import MPI
import h5py
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm
from numpy import fft
import matplotlib.ticker as mtick
import matplotlib
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from astropy.cosmology import LambdaCDM
from Fourier_Quad import Fourier_Quad



M = 1
m = 2
v = numpy.pi
l = 4
w = v/l
t = numpy.linspace(0,16,201)
M_x = M*v/(M+m)*t + m*l/(M+m)*numpy.sin(w*t)
M_y = -m*l/(M+m)*numpy.cos(w*t)

m_x = M*v/(M+m)*t - M*l/(M+m)*numpy.sin(w*t)
m_y = M*l/(M+m)*numpy.cos(w*t)

img = plot_tool.Image_Plot()
img.subplots(1,1)
img.axs[0][0].plot(m_x,m_y,marker="o",label="m")
img.axs[0][0].plot(M_x,M_y,marker="o",label="M")
# img.axs[0][1].scatter(t,numpy.sqrt((M_x-m_x)**2 + (M_y-m_y)**2))
img.axs[0][0].legend()
img.save_img("E:/Mm.png")
img.show_img()
exit()


ch_num = 8
cuts_num = 20
x_coord = [i*2/cuts_num*100 for i in range(ch_num)]
print(x_coord)
ch = [i*2 for i in range(ch_num)]
ylabels = ["m$_1 \\times 10^2$", "m$_2 \\times 10^2$", "m$_1 \\times 10^2$", "m$_2 \\times 10^2$"]
fmt = '%2.f%%'
xticks = mtick.FormatStrFormatter(fmt)

npz = numpy.load("E:/works/CFHT_tomo/all/cut_ext/flux_alt_s12_a1/total.npz")
mc1 = npz["arr_0"][:,ch]
mc2 = npz["arr_1"][:,ch]

img = plot_tool.Image_Plot()
img.subplots(1, 2)
img.axs[0][0].errorbar(x_coord, mc1[0]-1, mc1[1],marker="s", mfc="none",linewidth=img.plt_line_width,capsize=img.cap_size, label="$m_1$")
img.axs[0][0].errorbar(x_coord, mc2[0]-1, mc2[1], marker="s", mfc="none",linewidth=img.plt_line_width,capsize=img.cap_size, label="$m_2$")


img.axs[0][1].errorbar(x_coord, mc1[2], mc1[3], marker="s", mfc="none",linewidth=img.plt_line_width,capsize=img.cap_size, label="$c_1$")
img.axs[0][1].errorbar(x_coord, mc2[2], mc2[3], marker="s", mfc="none",linewidth=img.plt_line_width,capsize=img.cap_size, label="$c_2$")
for i in range(2):
    img.axs[0][i].xaxis.set_major_formatter(xticks)
    xs = img.axs[0][i].set_xlim()
    img.axs[0][i].plot([xs[0], xs[1]],[0,0],linestyle="--",c="grey")
    img.axs[0][i].legend(fontsize=img.legend_size)
img.set_label(0,0,0,"m")
img.set_label(0,1,0,"c")
img.set_label(0,0,1,"Cutoff percentage")
img.set_label(0,1,1,"Cutoff percentage")
img.subimg_adjust(h=0,w=0.26)
# img.axs[0][1].legend(fontsize=img.legend_size)
img.save_img("E:/cfht_cut.png")
img.show_img()
# seed = numpy.random.randint(1, 10000, 100)
# pool = []
# for i in range(10):
#     for j in range(10):
#         fq = Fourier_Quad(60,seed[i*10+j])
#         pts_num = 300
#         rand_pts = fq.ran_pts(num=pts_num, radius=9, ellip=0.8)
#         img_s = fq.convolve_psf(rand_pts, 4, 100, "Moffat")
#         pool.append(img_s)
# gal = fq.stack(pool,10)
my,mx = numpy.mgrid[0:48, 0:48]
fq = Fourier_Quad(48,123)
noise = fq.draw_noise(0,60)
gal = tool_box.gauss_profile(48,8,24,24,0.8,0.7)*200000 + noise
gal = tool_box.gauss_profile(48,3,24,24)*30000 + gal
gal_p = numpy.sqrt(fq.pow_spec(gal))
# img = plot_tool.Image_Plot()
# img.subplots(1, 2)
fig = plt.figure(figsize=plt.figaspect(0.5))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(mx, my, gal, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax = fig.add_subplot(1, 2, 2, projection='3d')
ax.plot_surface(mx, my, gal_p, rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
plt.show()
# img.show_img()
# img.close_img()
# print(seed.reshape((10,10)))


exit()

h5f = h5py.File("E:/ggl_test_result.hdf5","r")
print(list(h5f.keys()))
data = h5f["/data_0"].value
num = 5691
gt = data[0,:num]
gx = data[1,:num]
n = data[2,:num]
u = data[3,:num]
crit = data[4,:num]
g1 = data[5,:num]

gt = gt*crit
nt = (n+u)

print(gt)
print(g1)
print(crit)
# plt.hist(crit, 50)
# plt.hist(g1, 50)
plt.hist(gt, 1000)
plt.show()



exit()

cos = LambdaCDM(H0=70, Om0=0.31, Ode0=0.69)
r = cos.comoving_distance(0.4743854)
r1 = 2.99792458*1e5/70*0.4198876
print(r*0.7,r1*0.7)
exit()
def psf(flux, psf_scale, size, ellip, theta):
    my, mx = numpy.mgrid[0:size, 0:size] - size/2.
    r_scale_sq = 9
    m = 3.5

    rot_1 = numpy.cos(theta)
    rot_2 = numpy.sin(theta)
    q = (1+ellip)/(1-ellip)

    # mx_r = mx*r1 + my*r2
    # my_r = -mx*r2 + my*r1

    # factor = flux * 1. / (numpy.pi * psf_scale ** 2 * ((1. + r_scale_sq) ** (1. - m) - 1.) / (1. - m))
    # rsq = (mx_r/ psf_scale) ** 2 + (my_r/ psf_scale/q) ** 2
    # idx = rsq > r_scale_sq
    # rsq[idx] = 0.
    # arr = factor * (1. + rsq) ** (-m)
    # arr[idx] = 0.

    cent = size/2
    arr = numpy.zeros((size, size))
    rd = 1./psf_scale/psf_scale
    for i in range(size):
        ry1 = rot_2*(i-cent)
        ry2 = rot_1*(i-cent)
        for j in range(size):
            r1 = rot_1 * (j - cent) + ry1
            r2 = - rot_2 * (j - cent) + ry2
            rs = r1 * r1 * rd + r2 * r2 * rd*q

            if rs <= 9.:
                arr[i,j] += (1. + rs)**(-3.5)
                # arr[i, j] += numpy.exp(-rs*0.5)
    return arr

plt.subplot(121)
psf_img = psf(1, 30, 200, 0.7, 0)
idx = psf_img > 0
psf_img[idx] = 1
plt.imshow(psf_img)
plt.subplot(122)
psf_img = psf(1, 30, 200, 0.7, numpy.pi/4)
idx = psf_img > 0
psf_img[idx] = 1
plt.imshow(psf_img)
plt.show()


exit()
def plot_examples(cms):
    """
    helper function to plot two colormaps
    """
    numpy.random.seed(19680801)
    data = numpy.random.randn(30, 30)

    fig, axs = plt.subplots(1, 2, figsize=(6, 3), constrained_layout=True)
    for [ax, cmap] in zip(axs, cms):
        psm = ax.pcolormesh(data, cmap=cmap, rasterized=True, vmin=-4, vmax=4)
        fig.colorbar(psm, ax=ax)
    plt.show()

N = 256
vals = numpy.ones((N, 4))
vals[:, 0] = 0.95
vals[:, 1] = numpy.linspace(0, 0.01, N)
vals[:, 2] = numpy.linspace(0, 0.04, N)
vals[:, 3] = numpy.linspace(0, 1, N)
print(vals)
newcmp = ListedColormap(vals)

# Create custom colormaps
cdict = {'red': ((0.0, 1.0, 1.0),   # Full red at the first stop
                 (0.5, 0.0, 0.0),   # No red at second stop
                 (1.0, 1.0, 1.0)),  # Full red at final stop
        #
        'green': ((0.0, 0.0, 0.0),  # No green at all stop
                 (0.5, 0.0, 0.0),   #
                 (1.0, 0.0, 0.0)),  #
        #
        'blue': ((0.0, 0.0, 0.0),   # No blue at first stop
                 (0.5, 1.0, 1.0),   # Full blue at second stop
                 (1.0, 0.0, 0.0))}  # No blue at final stop

# cmap = LinearSegmentedColormap('Rd_Bl_Rd', cdict, 256)
# im = numpy.outer(numpy.ones(10), numpy.linspace(0, 255, 256))
# fig = plt.figure(figsize=(9, 2))
# ax = fig.add_subplot('111')
# ax.set_xticks(numpy.linspace(0, 255, 3))
# ax.set_xticklabels([0, 0.5, 1])
# ax.set_yticks([])
# # ax.set_yticklabels([])
im = numpy.abs(numpy.random.normal(0,100,10000).reshape((100,100)))
# norm = matplotlib.colors.Normalize(vmin=im.min(), vmax=im.max())
# ax_=ax.imshow(im, interpolation='nearest', cmap=newcmp)
#
# plt.colorbar(ax_, cmap=newcmp, norm=norm)
# plt.show()
img = plot_tool.Image_Plot(fig_x=8, fig_y=6)
img.set_style()
img.plot_img(1,1)
img.imgshow(0,0,im,True)
plt.show()
exit(0)

nx =1
ny =1
sub_fig = [[] for i in range(ny)]
print(sub_fig)
fig = plt.figure()
for i in range(ny):
    for j in range(nx):
        ax = fig.add_subplot(ny, nx, i * nx + j + 1)
        sub_fig[i].append(ax)
        print(type(sub_fig[i][j]))
        sub_fig[i][j].tick_params(direction='in')
print(sub_fig)
print(sub_fig[0][0])
exit()

f = h5py.File("E:/total.hdf5")
data = f["/mc1"].value
print(data)
print(list(f.keys()))


exit(0)
def pow_spec(image):
    image_ps = fft.fft2(image)
    return image_ps

def inv_pow(image):
    image_ps = fft.ifft2(image)
    return image_ps

def err(message):
    err_message = "%s is wrong"%message
    raise ValueError(err_message)

def mag(snr, r):
    return -numpy.log10(snr*numpy.sqrt(numpy.pi*(r/0.187)**2)*120)/0.4 + 34.5358
r = numpy.linspace(0.2, 1.2, 100)
plt.plot(r, mag(5, r))
plt.show()
exit()
x, y = numpy.mgrid[24-15:24+15,24-15:24+15]

fq = Fourier_Quad(48, 123)
gals = fq.segment(img)
# gal = gals[115]
my, mx = numpy.mgrid[0:48,0:48]
gal = numpy.exp(-((my-24)**2+(mx-24)**2)/2./25.)*10#+\
# hdu = fits.PrimaryHDU(gal)
# hdu.writeto("E:/gal.fits",overwrite=True)
     #numpy.exp(-((my-20)**2+(mx-24)**2)/2/100)*100
# points = fq.ran_pos(num=100, radius=20, g=(0.01, 0))[1]
# gal = fq.convolve_psf(pos=points, psf_scale=4, flux=1000, psf="Moffat")# + noise


gal_p = numpy.log10(fq.pow_spec(gal))
# hdu = numpy.log10(gal_p)
# smooth_hdu = tool_box.smooth(hdu,48)
# # w = fits.PrimaryHDU(hdu)
# # w.writeto("E:/test.fits",overwrite=True)
# test = fits.open("E:/smooth.fits")[0].data
# plt.subplot(131)
# plt.imshow((test-smooth_hdu)/test)
# plt.colorbar()
# plt.subplot(132)
# plt.imshow(smooth_hdu)
# plt.colorbar()
# plt.subplot(133)
# plt.imshow(test)
# plt.colorbar()
# plt.show()
# exit()
gal_ps = tool_box.smooth(gal_p, 48)

peak_pow = numpy.max(gal_p)
fit_area = numpy.ones((5,5))*peak_pow
fit_x, fit_y = numpy.mgrid[22:27, 22:27]
rim = fq.border(1)
num = numpy.sum(rim)
fsnr = numpy.sqrt(10**gal_p[24,24]*num/numpy.sum(10**gal_p*rim))
fsnr_f = numpy.sqrt(10**gal_ps[24,24]*num/numpy.sum(10**gal_p*rim))
print(fit_area.shape, fit_x.shape, gal_p.shape, x.shape)
print(gal_ps.max()/gal_p.max(), fsnr, fsnr_f)

fig = plt.figure(figsize=(14,14))
ax1 = fig.add_subplot(231)
ax1.imshow(gal)
ax2 = fig.add_subplot(232)
ax2.set_title("Powerspectrum, SNR=%.2f"%fsnr)
ax2.imshow((gal_p - gal_ps)/gal_p)
ax3 = fig.add_subplot(233)
ax3.imshow(gal_ps)
ax3.set_title("Smoothed powerspectrum, SNR$_{fit}$=%.2f"%fsnr_f)
ax4 = fig.add_subplot(234)
ax4.imshow(gal_p)
# surf = ax4.plot_surface(x, y, gal_p[24-15:24+15,24-15:24+15], rstride=1, cstride=1, cmap=cm.coolwarm,
#                        linewidth=0, antialiased=False)
# fig.colorbar(surf, shrink=0.5, aspect=10)
# ax4.plot_surface(fit_x, fit_y, fit_area)

ax5 = fig.add_subplot(235, projection='3d')
surf = ax5.plot_surface(x, y, gal_p[24-15:24+15,24-15:24+15], rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=0.5, aspect=10)
ax5.plot_surface(fit_x, fit_y, fit_area)
ax5.set_title("Peak: %g"%peak_pow)
ax6 = fig.add_subplot(236, projection="3d")
ax6.plot_surface(x, y, gal_ps[24-15:24+15,24-15:24+15], rstride=1, cstride=1, cmap=cm.coolwarm,
                       linewidth=0, antialiased=False)
ax6.set_title("Peak: %g"%gal_ps.max())
plt.suptitle("$\\frac{SNR}{SNR_{fit}}$ = %g"%(gal_p.max()/gal_ps.max()), fontsize=20)
plt.subplots_adjust(wspace=0.1,hspace=0.1)
plt.savefig("E:/2.png",bbox_inches='tight')
plt.show()


# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# cpus = comm.Get_size()
#
# data_path = "/mnt/ddnfs/data_users/hkli/simu_test_gal_bigger/result/data/data_2.0sig/data_%d.hdf5"%rank
# f = h5py.File(data_path)
# data = f["/data"].value
# f.close()
# snr_f = data[:, 2]
# snr = data[:, 3]
# mag = data[:, 6]
# idx1 = mag< 25
# idx2 = mag>24.5
# idx = data[:, 1] > 0
# plt.figure(figsize=(12, 6))
# plt.subplot(221)
# plt.scatter(mag[idx], snr[idx], s=0.3)
# plt.ylim(10**(-4),10**2.5)
# plt.yscale("log")
# plt.subplot(222)
# plt.scatter(mag[idx], snr_f[idx], s=0.3)
# plt.ylim(10**(-4),10**2.5)
# plt.yscale("log")
# plt.subplot(223)
# plt.scatter(snr[idx], snr_f[idx], s=0.3)
# # plt.hist(x=snr[idx],bins=40)
# plt.subplot(224)
# plt.hist(x=snr_f[idx], bins=40)
#
# plt.savefig("/home/hkli/work/test/snr_hist_%d.png"%rank)
# plt.close()
#




