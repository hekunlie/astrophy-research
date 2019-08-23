from sys import path
path.append("E:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import numpy
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.axes3d import Axes3D, get_test_data
from matplotlib import cm



def plot_3d(fig, row, col, seq, x, y, z, labels):
    label_size = 20
    ax = fig.add_subplot(row, col, seq, projection='3d')
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    ax.set_xlabel(labels[0], fontsize=label_size)
    ax.set_ylabel(labels[1], fontsize=label_size)
    ax.set_zlabel(labels[2], fontsize=label_size)
    fig.colorbar(surf, shrink=0.5, aspect=10)

def find_patch(arr,thresh):
    ny, nx = arr.shape
    x1, x2 = -1, -1
    y1, y2 = -1, -1
    for i in range(int(ny/2)):
        if arr[i].min() <= thresh and y1 < 0:
            y1 = i
        if arr[ny-i-1].min() <= thresh and y2 < 0:
            y2 = ny-i
    for i in range(int(nx/2)):
        if arr[:,i].min() <= thresh and x1 < 0:
            x1 = i
        if arr[:, nx-i-1].min() <= thresh and x2 < 0:
            x2 = nx-i
    print(y1, y2, x1, x2)

    if min(x1, x2, y1, y2) >= 0 and x2-x1 > 1 and y2-y1>1:
        return y1, y2, x1, x2
    else:
        return 0, ny, 0, nx



gal_num = 100000
x = numpy.random.uniform(-8, 8, gal_num)
a1, a2 = 0, -0.01
shear_slope_1d = a1 + a2*x
print(shear_slope_1d.max())
fq = Fourier_Quad(6,12)

bin_num = 8
bin_num2 = int(bin_num/2)
inverse = range(int(bin_num2 - 1), -1, -1)

bin_num_test = 8
bin_num2_test = int(bin_num_test/2)
inverse_test = range(int(bin_num2_test - 1), -1, -1)

hist_bin = 20

# simulation
mg = numpy.zeros((gal_num,))
mnu = numpy.ones_like(mg)
for i in range(gal_num):
    mg[i] = numpy.random.normal(shear_slope_1d[i], 0.3, 1)
mg_bins = fq.set_bin(mg, bin_num, 10)
mg_bins_test = fq.set_bin(mg, bin_num_test, 10)

num_ini = numpy.histogram(mg, mg_bins)[0]
n1 = num_ini[0:bin_num2][inverse]
n2 = num_ini[bin_num2:]
num_exp = (n1 + n2) / 2


img = Image_Plot(fig_x=12, fig_y=9)
img.subplots(1, 1)
chi_temp = fq.get_chisq(mg, mnu, 1, mg_bins, bin_num2, inverse, 0)
img.axs[0][0].hist(mg, hist_bin, histtype='step',linewidth=img.plt_line_width+2, label="Before PDF, $\chi^2$=%.2f"%chi_temp)

# num_null = numpy.histogram(mg, mg_bins_test)[0]
# n1 = num_null[0:bin_num2_test][inverse_test]
# n2 = num_null[bin_num2_test:]
# chi_null = (n1-n2)/numpy.abs(n1-n2) * numpy.sum((n1-n2)**2/(n1+n2))*0.5
# img.axs[0][1].scatter(range(bin_num2_test), chi_null,marker="p", s=140, label="Before PDF, $\chi^2$=%.2f"%chi_temp)

a2_test = numpy.linspace(-0.01-0.1, -0.01+0.1, 5)
print(a2_test)
for i in range(5):
    gh = - a2_test[i]*x
    chi_temp = fq.get_chisq(mg, mnu, gh, mg_bins, bin_num2, inverse, 0)
    img.axs[0][0].hist(mg-gh, hist_bin, histtype='step',linewidth=img.plt_line_width+2, label="$a_2$=%.2f, $\chi^2$=%.2f"%(a2_test[i],chi_temp))

    # num_test = numpy.histogram(mg-gh, mg_bins_test)[0]
    # n1 = num_test[0:bin_num2_test][inverse_test]
    # n2 = num_test[bin_num2_test:]
    # chi_test = (n1-n2)/numpy.abs(n1-n2) * numpy.sum((n1-n2)**2/(n1+n2))*0.5
    # img.axs[0][1].scatter(range(bin_num2_test), chi_test,s=120, label="$a_2$=%.2f, $\chi^2$=%.2f"%(a2_test[i],chi_temp))
    # img.axs[0][2].plot(range(bin_num2_test), chi_test,linewidth=img.plt_line_width, label="$\hat g$=%.2f, $\chi^2$=%f"%(a2_test[i],chi_temp))
ys = img.axs[0][0].set_ylim()
img.axs[0][0].plot([0,0],[ys[0], ys[1]], linestyle="--", c="grey")
img.axs[0][0].set_ylim(ys)
img.axs[0][0].legend(fontsize=img.legend_size-2)
img.show_img()



exit()

nx, ny = 100, 100

a1_range = numpy.linspace(-0.1, 0.1, nx)
a2_range = numpy.linspace(-0.1, 0.1, ny)

a1_hat, a2_hat = numpy.meshgrid(a1_range, a2_range)
chisq = numpy.zeros((ny, nx))
chisq_new = numpy.zeros((ny, nx))

# for i in range(ny):
#     for j in range(nx):
#         gh = a1_hat[i,j] + a2_hat[i, j]*x
#         chisq[i,j] = fq.get_chisq(mg, mnu, gh, mg_bins, bin_num2,inverse, 0)
#         chisq_new[i,j] = fq.get_chisq_new(mg, mnu, gh, mg_bins, bin_num2, inverse, 0, num_exp)
# numpy.savez("./chisq.npz", chisq, chisq_new)
f = numpy.load("./chisq.npz")
chisq = f["arr_0"]
chisq_new = f["arr_1"]

y1,y2,x1, x2 = find_patch(chisq, 0)
y1_n,y2_n,x1_n, x2_n = find_patch(chisq_new, 0)

# img = Image_Plot()
# img.subplots(1,2)
# sub_fig = img.axs[0][0].imshow(chisq[y1:y2,x1:x2])
# plt.colorbar(sub_fig, ax=img.axs[0][0],cmap="jet")
# sub_fig = img.axs[0][1].imshow(chisq_new[y1_n:y2_n,x1_n:x2_n])
# plt.colorbar(sub_fig, ax=img.axs[0][1],cmap="jet")
# img.show_img()

# fig = plt.figure(figsize=(16, 8))
# plot_3d(fig, 1,2,1,a1_hat, a2_hat, chisq,["$a_1$","$a_2$","$\chi^2$"])
# plot_3d(fig, 1,2,2,a1_hat, a2_hat, chisq_new,["$a_1$","$a_2$","$\chi^2_{new}$"])
# plot_3d(fig, 1,2,1,a1_hat[y1:y2,x1:x2], a2_hat[y1:y2,x1:x2], chisq[y1:y2,x1:x2],["$a_1$","$a_2$","$\chi^2$"])
# plot_3d(fig, 1,2,2,a1_hat[y1_n:y2_n,x1_n:x2_n], a2_hat[y1_n:y2_n,x1_n:x2_n], chisq_new[y1_n:y2_n,x1_n:x2_n],
#         ["$a_1$","$a_2$","$\chi^2_{new}$"])
# plt.show()