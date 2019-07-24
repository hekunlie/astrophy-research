import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
import h5py
import MCMC_program



parent_path="/mnt/perc/hklee/CFHT/multi_shear/cluster_field/"

fq = Fourier_Quad(10,123)
for i in range(50):
    h5f = h5py.File(parent_path + "result/expo_%d_slope.hdf5"%i,"r")
    temp = h5f["/data"].value
    h5f.close()
    if i == 0:
        data = temp
    else:
        data = numpy.row_stack((data, temp))

mg1 = data[:, 0]
mg2 = data[:, 1]
mnu1 = data[:, 2] + data[:, 3]
mnu2 = data[:, 2] - data[:, 3]
ra = data[:, 6]
dec = data[:, 7]

nx, ny = 20, 20
delta_ra = 16 #arcmin
delta_dec = delta_ra
half_side = nx/2
# arcmin/pix
pixel_scale = delta_ra/nx

ra_bin = numpy.linspace(-delta_ra/2, delta_ra/2, nx+1)
dec_bin = numpy.linspace(-delta_dec/2, delta_dec/2, ny+1)

shear_field_cal = numpy.zeros((nx, ny))
num_in_grid = numpy.zeros((nx, ny))
for i in range(ny):
    idx1 = dec >= dec_bin[i]
    idx2 = dec <= dec_bin[i + 1]
    for j in range(nx):
        idx3 = ra >= ra_bin[j]
        idx4 = ra <= ra_bin[j+1]

        idx = idx1 & idx2 & idx3 & idx4

        num_in_grid[i,j] = idx.sum()
        gh, gh_sig = fq.fmin_g_new(mg1[idx], mnu1[idx], 8)[:2]
        shear_field_cal[i, j] = gh
        print("%.4f, %.4f"%(gh, gh_sig))



# the grid
# x: ra, y: dec
my, mx = numpy.mgrid[-half_side:half_side,-half_side:half_side]
y, x = my*pixel_scale, mx*pixel_scale

ra_bin = numpy.linspace(-delta_ra/2, delta_ra/2, nx+1)
dec_bin = numpy.linspace(-delta_dec/2, delta_dec/2, ny+1)

ra_range = (-delta_ra/2, delta_ra/2)
dec_range = (-delta_dec/2, delta_dec/2)

a1, a2, a3 = 0, 0.002, -0.006
profile_params = [a1, a2, a3]
# shear field
shear_field = MCMC_program.shear_slope(profile_params, x, y)

img = Image_Plot()
img.subplots(1, 3)
img.set_style()
inverse = range(nx - 1, -1, -1)
show_data = [shear_field[inverse], shear_field_cal[inverse], num_in_grid[inverse]]
titles = ["Input", "Recovered", "Num"]

g_max = max(show_data[0].max(), show_data[1].max())
g_min = max(show_data[0].min(), show_data[1].min())
norm_g = plt.Normalize(vmin=g_min, vmax=g_max)
cmap_g = plt.get_cmap('jet')

figs = []
for i in range(3):
    # if i < 2:
    #     fig = img.axs[0][i].imshow(show_data[i], cmap="jet")
    #     # # sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
    #     # # sm._A = []
    #     # plt.colorbar(fig, ax=img.axs[0][i])
    # else:
    #     fig = img.axs[0][i].imshow(show_data[i], cmap="jet")
    #     img.figure.colorbar(fig, ax=img.axs[0][i])

    # figs.append(fig)

    fig = img.axs[0][i].imshow(show_data[i], cmap="jet")
    img.figure.colorbar(fig, ax=img.axs[0][i])

    img.set_label(0, i, 0, "DEC. [arcmin]")
    img.set_label(0, i, 1, "R.A. [arcmin]")
    if i < 2:
        pass

    fmt = ["%.1f"%(8 - 4*k) for k in range(5)]
    img.set_ticklabel(0, i, 0, len(fmt), fmt)
    fmt = ["%.1f"%(-8 + 4*k) for k in range(5)]
    img.set_ticklabel(0, i, 1, len(fmt), fmt)

    # img.axs[0][i].set_yticks([int(nx - 1 - k * nx / 5) for k in range(6)])
    # img.axs[0][i].set_yticklabels(["%.1f" % (k * nx / 5. * pixel_scale) for k in range(6)])
    # img.axs[0][i].set_xticks([int(k * ny / 5) for k in range(6)])
    # img.axs[0][i].set_xticklabels(["%.1f" % (k * ny / 5. * pixel_scale) for k in range(6)])

    img.axs[0][i].set_title(titles[i])

img.save_img(parent_path + "pic/shear_field_slope_test.png")
img.show_img()


#

#
# gh, gh_sig = fq.fmin_g_new(mg1, mnu1, 8)[:2]
#
# print(gh, gh_sig)
# gh, gh_sig = fq.fmin_g_new(mg2, mnu2, 8)[:2]
#
# print(gh, gh_sig)

# for i in range(30):
#     h5f = h5py.File(parent_path + "result/expo_%d_slope.hdf5"%i,"r")
#     data = h5f["/data"].value
#     h5f.close()
#
#     mg1 = data[:, 0]
#     mg2 = data[:, 1]
#     mnu1 = data[:, 2] + data[:, 3]
#     mnu2 = data[:, 2] - data[:, 3]
#
#     gh, gh_sig = fq.fmin_g_new(mg1, mnu1, 8)[:2]
#
#     print("%.5f, %.5f, %.5f"%(0.002*i, gh, gh_sig))