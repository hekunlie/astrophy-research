import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import tool_box
import correlation_function_tool as cf_tool
import numpy
from plot_tool import Image_Plot
import matplotlib.pyplot as plt
import time
import scipy
import h5py


# parameters
c = 299792.458  # Km/s
H0 = 67.5 #Km/s/Mpc
h = H0 / 100  # Km/s/Mpc
As = 2e-9
ns = 0.965

omg_cm0 = 0.26
omg_cm0h2 = 0.26*h**2

omg_bm0 = 0.05
omg_bm0h2 = 0.05*h**2

omg_m0 = omg_cm0+omg_bm0

alpha = 9 / 4 * omg_m0 ** 2 * (100/c) ** 4  # h^4 \cdot Mpc^{-4}
print(alpha)

kmin, cal_kmax, interp_kmax, kpts_num = 1e-4, 3, 4, 300

zmin, interp_zmax = 0, 4

Lpts_min, Lpts_max, Lpts_num = 10, 8000, 10000

theta_num  = 50
theta_min, theta_max = 0.8, 60 # arcmin
theta = numpy.linspace(theta_min, theta_max,theta_num)
theta_radian = theta/60/180*numpy.pi


redshift_bin_num = 6
redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3],dtype=numpy.float32)

nz_bin_num = 400

# read data
h5f = h5py.File("/mnt/perc/hklee/CFHT/correlation/cata/stack_data.hdf5","r")
data = h5f["/data"][()]
h5f.close()


ts = time.time()
# tomo panel num
tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)

# histogram
zehist, zebin, zebin_cent = cf_tool.get_nz(data[:, 8], redshift_bin, data[:, 9], nz_bin_num)

inv_scale_factor_sq = (1 + zebin_cent) ** 2

# comoving distance Mpc/h
com_dist = cf_tool.z2dist(zebin_cent, H0, omg_m0)
delta_com_dist = com_dist[1:] - com_dist[:-1]

# lense efficiency
qx = cf_tool.get_lenq(zehist, com_dist)

# g(x)^2/a(x)^2 in the integration
integ_factor = numpy.zeros((tomo_panel_num, nz_bin_num))
tag = 0
for i in range(redshift_bin_num):
    for j in range(i, redshift_bin_num):
        integ_factor[tag] = qx[i] * qx[j] * inv_scale_factor_sq
        tag += 1

# set up bins for L of P(L), and decide the ks needed
Lpts, dLpts, integ_kh = cf_tool.ready4PL(Lpts_min, Lpts_max, Lpts_num, kmin, cal_kmax, com_dist, nz_bin_num)

Lpts_theta = numpy.zeros((theta_num, Lpts_num))

for i in range(theta_num):
    Lpts_theta[i] = theta_radian[i] * Lpts

# Pk interpolation
integ_pk = numpy.zeros((nz_bin_num, Lpts_num))
PLs = numpy.zeros((tomo_panel_num, Lpts_num))
# the theoretical line
xi_plus = numpy.zeros((tomo_panel_num, theta_num))

zfit = numpy.linspace(zmin, interp_zmax, 150)

camb_result = cf_tool.get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, zfit, kmax=interp_kmax)[0]
pk_interp = camb_result.get_matter_power_interpolator(nonlinear=True, hubble_units=True, k_hunit=True)
sigma8_now = camb_result.get_sigma8()[-1]
print("sigma8: %.4f" % sigma8_now)

for i in range(nz_bin_num):
    idx = integ_kh[i] > 0
    integ_pk[i][idx] = pk_interp.P(zebin_cent[i], integ_kh[i][idx])

# get P(L)
for i in range(tomo_panel_num):
    pls_i = cf_tool.get_PL(integ_pk, Lpts_num, integ_factor[i], delta_com_dist) * alpha
    PLs[i] = pls_i

    # calculate the \chi_plus(\theta) , the theoretical line
    for j in range(theta_num):
        integ_part = Lpts * scipy.special.j0(Lpts_theta[j]) * PLs[i]
        xi_plus[i, j] = numpy.sum(((integ_part[1:] + integ_part[:-1]) / 2 * dLpts)) / 2 / numpy.pi

te = time.time()
print(te - ts)

numpy.savez("tomo_pk.npz", qx, integ_factor, PLs, xi_plus)


labels = ["%d-%d" % (i, j) for i in range(1, redshift_bin_num + 1) for j in range(i, redshift_bin_num + 1)]

img = Image_Plot(xpad=0.25, ypad=0.25)
img.subplots(2, 3)

for i in range(redshift_bin_num):
    img.axs[0][0].plot(zebin_cent, zehist[i], label="z bin %d" % i)
    img.axs[0][1].plot(com_dist, qx[i], label="z bin %d" % i)

for i in range(tomo_panel_num):
    ls, c = img.line_styles[i]
    img.axs[0][2].plot(com_dist, integ_factor[i], ls=ls, c=c, label=labels[i])

img.set_label(0, 1, 0, "g($\chi$)")
img.set_label(0, 1, 1, "$\chi [Mpc\cdot h^{-1}]$")

img.set_label(0, 2, 0, "$g_i(\chi)g_j(\chi)/a(\chi)^2$")
img.set_label(0, 2, 1, "$\chi [Mpc\cdot h^{-1}]$")
for i in range(3):
    img.axs[0][i].legend(ncol=3)

for i in range(tomo_panel_num):
    ls, c = img.line_styles[i]
    img.axs[1][0].plot(Lpts, Lpts * (Lpts + 1) * PLs[i] / numpy.pi / 2, label=labels[i], ls=ls, c=c)
    img.axs[1][1].plot(theta, chi_plus[i], label=labels[i], ls=ls, c=c)

# img.axs[0][0].set_xlim(0.7,10000)
img.axs[1][0].set_ylim(1e-10, 1e-3)
for i in range(2):
    img.axs[1][i].set_xscale("log")
    img.axs[1][i].set_yscale("log")

img.set_label(1, 0, 0, "$\ell(\ell+1)}$/2$\pi$ $P_{\kappa}(\ell)$")
img.set_label(1, 0, 1, "$\ell$")

img.set_label(1, 1, 0, "$\\xi_{+}(\\theta)$")
img.set_label(1, 1, 1, "$\\theta$")

img.axs[1][0].legend(ncol=3)
img.axs[1][1].legend(ncol=3)

img.figure.delaxes(img.axs[1][2])
img.save_img("tomo_pk.png")