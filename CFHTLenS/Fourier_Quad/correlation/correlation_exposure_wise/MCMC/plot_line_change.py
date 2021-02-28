from sys import path,argv
path.append("/home/hklee/work/mylib")
import time
from plot_tool import Image_Plot
import tool_box
import numpy
import corner
import pyccl as ccl
import correlation_function_tool as cf_tool
import h5py
import camb
import time

cmd = argv[1]
dx = float(argv[2])
num = int(argv[3])

Lpts_min, Lpts_max, Lpts_num = 10, 20000, 10000
ell = tool_box.set_bin_log(Lpts_min, Lpts_max, Lpts_num)

redshift_bin = numpy.array([0.2, 0.39, 0.58, 0.72, 0.86, 1.02, 1.3], dtype=numpy.float32)
redshift_bin_num = len(redshift_bin) - 1
tomo_panel_num = int((redshift_bin_num * redshift_bin_num + redshift_bin_num) / 2)

theta_num = 30
theta = numpy.zeros((tomo_panel_num, theta_num))
theta_radian = numpy.zeros((tomo_panel_num, theta_num))

for i in range(tomo_panel_num):
    theta[i] = tool_box.set_bin_log(1, 60, theta_num)
    theta_radian[i] = theta[i] / 60 / 180 * numpy.pi
theta_deg = theta / 60

h5f = h5py.File("/mnt/perc/hklee/CFHT/correlation/cata/stack_data.hdf5", "r")
data = h5f["/data"][()]
h5f.close()

nz_bin_num = 340

zehist, zebin, zebin_cent = cf_tool.get_nz(data[:, 8], redshift_bin, data[:, 9], nz_bin_num, (0, 3), False)
# img = Image_Plot(xpad=0.25)
# img.subplots(1, 1)
# for j in range(redshift_bin_num):
#     img.axs[0][0].plot(zebin_cent, zehist[j], label="[%.2f~%.2f]" % (redshift_bin[j], redshift_bin[j + 1]))
# img.axs[0][0].legend()
# img.show_img()
# img.close_img()

c = 299792.458
H0 = 67.4
h = H0 / 100
ns = 0.965
sigma8_planck = 0.811
omega_c_planck = 0.2642
omega_b_planck = 0.0493

if cmd == "sig8":
    para_list = numpy.linspace(sigma8_planck - dx, sigma8_planck + dx, num)
elif cmd == "omega_c":
    para_list = numpy.linspace(omega_c_planck - dx, omega_c_planck + dx, num)
elif cmd == "omega_b":
    para_list = numpy.linspace(omega_b_planck - dx, omega_b_planck + dx, num)
else:
    print("Wrong para")
    exit()

img_pk = cf_tool.pre_img(redshift_bin_num)

for i in range(num):
    if cmd == "sig8":
        sigma8 = para_list[i]
        omega_c = omega_c_planck
        omega_b = omega_b_planck
    elif cmd == "omega_c":
        sigma8 = sigma8_planck
        omega_c = para_list[i]
        omega_b = omega_b_planck
    elif cmd == "omega_b":
        sigma8 = sigma8_planck
        omega_c = omega_c_planck
        omega_b = para_list[i]


    omega_m = omega_b + omega_c
    omega_ch2 = omega_c * h * h
    omega_bh2 = omega_b * h * h
    omega_mh2 = omega_m * h * h
    print(sigma8, omega_ch2, omega_bh2)

    t1 = time.time()
    ccl_xip, ccl_xim, ccl_PL = cf_tool.get_tomo_xi_ccl(sigma8, omega_c, omega_b, h, zebin_cent, zehist, theta_deg, ell)

    t2 = time.time()

    print(t2 - t1)

    #     cache_file = tool_box.file_name("/mnt/perc/hklee/CFHT/correlation/test/lines/line.hdf5")
    cache_file = "/mnt/perc/hklee/CFHT/correlation/test/lines/line_%s_%d.hdf5" %(cmd, i)

    h5f = h5py.File(cache_file, "w")
    h5f["/xi_p"] = ccl_xip
    h5f["/xi_m"] = ccl_xim
    h5f["/paras"] = numpy.array([sigma8, omega_c, omega_b, h])
    h5f["/paras_name"] = "sigma8, omega_c, omega_b,h"
    h5f.close()

    lb = "$\sigma_8$=%.3f, $\Omega_m$=%.3f\n$\Omega_b$=%.3f,ccl $\\xi_{+}$" % (sigma8, omega_c, omega_b)
    cf_tool.plot_panel(img_pk, redshift_bin_num, theta, ccl_xip, label=lb, color="C%d" % i, xlim=(0.8, 70), ylim=None)
    cf_tool.plot_panel(img_pk, redshift_bin_num, theta, ccl_xim, label="ccl $\\xi_{-}$", color="C%d" % i, ls="--",
                       xlim=(0.8, 70), ylim=None)

img_pk.save_img("/mnt/perc/hklee/CFHT/correlation/test/lines/line_%s.png"%cmd)
img_pk.close_img()
# img_camb.show_img()
