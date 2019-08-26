from sys import path
path.append("D:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm





seed = 12114
total_num = 25000
uni_num = 15000

a1, a2 = 0, -0.01

fq = Fourier_Quad(6,12)
rng = numpy.random.RandomState(seed)

bin_num = 8
bin_num2 = int(bin_num/2)
inverse = range(int(bin_num2 - 1), -1, -1)

bin_num_test = 8
bin_num2_test = int(bin_num_test/2)
inverse_test = range(int(bin_num2_test - 1), -1, -1)

hist_bin = 20

# x = numpy.random.uniform(-8, 8, gal_num)
x_asy = rng.uniform(0, 8, int(total_num - uni_num))
x_uni = rng.uniform(-8, 8, uni_num)
x = numpy.zeros((total_num,))

x[:uni_num] = x_uni
x[uni_num:] = x_asy

shear_slope_1d = a1 + a2*x
print(shear_slope_1d.min(), shear_slope_1d.max())

# img = Image_Plot(fig_x=8,fig_y=8)
# img.subplots(1,1)
# nx, ny = 71, 71
# a1_range = numpy.linspace(-0.1, 0.1, nx)
# a2_range = numpy.linspace(-0.1, 0.1, ny)
# for i in range(nx):
#     for j in range(nx):
#         img.axs[0][0].scatter(a1_range[i], a2_range[j],marker="s",s=10)
# img.axs[0][0].set_xlim((-0.1156491701936468, 0.11153579987303396))
# img.axs[0][0].set_ylim((-0.11758202200966764, 0.11207960506016676))
# img.show_img()
# exit()
# simulation
mg_ori = numpy.zeros((total_num,))
mg = numpy.zeros((total_num,))
mnu = numpy.ones_like(mg)
for i in range(total_num):
    mg_ori[i] = rng.normal(0, 0.3, 1)
    # mg[i] = rng.normal(shear_slope_1d[i], 0.3, 1)
    mg[i] = mg_ori[i] + shear_slope_1d[i]

mg_bins = fq.set_bin(mg, bin_num, 10)
mg_bins_test = fq.set_bin(mg, bin_num_test, 10)

num_ini = numpy.histogram(mg, mg_bins)[0]
n1 = num_ini[0:bin_num2][inverse]
n2 = num_ini[bin_num2:]
num_exp = (n1 + n2) / 2

test_num = 15
a2_test = numpy.linspace(a2-0.1, a2+0.1, test_num)

img = Image_Plot()
img.subplots(1,2)
x_plt = numpy.linspace(-8, 8, 20)
fx = a2*x_plt
# img.axs[0][1].plot(x_plt, fx)

chi_temp = fq.get_chisq(mg_ori, 0, 0, mg_bins, bin_num2, inverse, 0)
chi_temp_sym = fq.get_chisq(mg, 1, a2*x, mg_bins, bin_num2, inverse, 0)
print(chi_temp,chi_temp_sym)
bin_idx = []
for i in range(bin_num):
    idx_1 = mg_ori >= mg_bins[i]
    idx_2 = mg_ori < mg_bins[i+1]
    idx = idx_1 & idx_2
    bin_idx.append(idx)
    print(idx.shape)
    img.axs[0][0].scatter(mg_ori[idx], x[idx], c="C%d"%i,s=7)
    img.axs[0][1].scatter(mg[idx], x[idx], c="C%d"%i,s=7)


for i in range(2):
    ys = img.axs[0][i].set_ylim()
    x_label = numpy.linspace(-0.5,0.5,9)
    for j in range(9):
        img.axs[0][i].plot([x_label[j],x_label[j]],[ys[0], ys[1]],linestyle="--", c="k", linewidth=2)
    img.axs[0][i].set_ylim(ys)
img.show_img()


for k in range(test_num):
    mg_i = mg - a2_test[k]*x
    chi_temp_sym = fq.get_chisq(mg_i, 0, 0, mg_bins, bin_num2, inverse, 0)
    img = Image_Plot()
    img.subplots(1,2)
    for ib in range(bin_num):
        img.axs[0][0].scatter(mg_i[bin_idx[ib]], x[bin_idx[ib]], c="C%d" %ib, s=7)

    ys = img.axs[0][0].set_ylim()
    x_label = numpy.linspace(-0.5, 0.5, 7)
    for j in range(7):
        img.axs[0][0].plot([x_label[j], x_label[j]], [ys[0], ys[1]], linestyle="--", c="k", linewidth=2)
    img.axs[0][0].set_ylim(ys)
    img.axs[0][0].set_title("$\chi^2$=%.2f, $a_2$=%.3f"%(chi_temp_sym,a2_test[k]),fontsize=img.xy_lb_size)
    img.set_label(0,0,0,"X")
    img.set_label(0,0,1,"G")

    fx = a2_test[k]*x_plt
    print(fx.min(), fx.max())
    img.axs[0][1].plot(fx,x_plt)
    img.set_label(0,1,1,"g")
    img.set_label(0,1,0,"X")
    img.axs[0][1].set_xlim(-0.9,0.9)

    img.save_img("H:/GGL/1D_slope/move/%d.png"%k)
    # img.show_img()
    img.close_img()
