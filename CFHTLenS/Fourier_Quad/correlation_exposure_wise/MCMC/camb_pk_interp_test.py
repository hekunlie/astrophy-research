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



c = 299792.458  # Km/s
H0 = 70.5 #Km/s/Mpc
h = H0 / 100  # Km/s/Mpc
As = 2.4e-9
ns = 0.965

omg_cm0 = 0.26
omg_cm0h2 = 0.26*h**2

omg_bm0 = 0.05
omg_bm0h2 = 0.05*h**2

omg_m0 = omg_cm0+omg_bm0

kmin, cal_kmax, interp_kmax, kpts_num = 1e-4, 3, 5, 300

zmin, interp_zmax, zbin_num, interp_znum = 0, 4, 150, int(argv[1])

zpts_test_1 = numpy.linspace(zmin, interp_zmax, zbin_num)

zpts_test_2 = numpy.linspace(zmin, interp_zmax, interp_znum)

kh_test = tool_box.set_bin_log(kmin, cal_kmax, kpts_num)


t1 = time.time()

camb_result = cf_tool.get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, zpts_test_1, kmax=interp_kmax)[0]

pk_interp = camb_result.get_matter_power_interpolator()
sigma8_1 = camb_result.get_sigma8()[-1]
pk_nonlin_interp_1 = cf_tool.get_interpPk(pk_interp, zpts_test_1, kh_test)


t2 = time.time()

camb_result = cf_tool.get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, zpts_test_2, kmax=interp_kmax)[0]

pk_interp = camb_result.get_matter_power_interpolator()
sigma8_2 = camb_result.get_sigma8()[-1]
pk_nonlin_interp_2 = cf_tool.get_interpPk(pk_interp, zpts_test_1, kh_test)


t3 = time.time()

camb_result = cf_tool.get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, zpts_test_1, kmax=cal_kmax)[0]
z_nonlin, kh_nonlin, pk_nonlin, s8 = cf_tool.get_PK(camb_result, kpts_num=kpts_num, kmax=cal_kmax)

t4 = time.time()

diff = kh_test - kh_nonlin

print("Time: %.3f  %.3f  %.3f"%(t2 - t1, t3 - t2, t4 - t3))
print("sigma8: %.5f  %.5f  %.5f"%(sigma8_1, sigma8_2, s8[-1]))
print(numpy.max(numpy.abs(diff)))

norm = plt.Normalize(vmin=numpy.min(zpts_test_1), vmax=numpy.max(zpts_test_1))
cmap = plt.get_cmap('YlOrRd')
cl = cmap(norm(zpts_test_1))

img = Image_Plot(fig_x=8, fig_y=5)
img.subplots(1, 2)

for tag, i in enumerate(zpts_test_1):
    pk_test = pk_nonlin_interp_1[tag, :]
    pk_nonlin_test = pk_nonlin[tag, :]
    diff = (pk_test - pk_nonlin_test) / pk_nonlin_test * 100

    img.axs[0][0].plot(kh_nonlin, diff, c=cl[tag])

    pk_test = pk_nonlin_interp_2[tag, :]
    pk_nonlin_test = pk_nonlin[tag, :]
    diff = (pk_test - pk_nonlin_test) / pk_nonlin_test * 100

    img.axs[0][1].plot(kh_nonlin, diff, c=cl[tag])

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []

img.axs[0][0].set_title("Pk Interpolation test, same Z")
img.axs[0][1].set_title("Pk Interpolation test, diff Z")

for i in range(2):
    clb = img.figure.colorbar(sm, ax=img.axs[0][i])
    clb.ax.set_xlabel("Z")

    img.set_label(0, i, 0, "(interp Pk - true Pk)/true Pk [%]")
    img.set_label(0, i, 1, "k [$ h\cdot Mpc^{-1}]$")
    img.axs[0][i].set_xscale("log")
    img.axs[0][i].set_xlim(7e-5, 4)
img.save_img("pk_interp_test_%dz.png"%interp_znum)
# img.show_img()
