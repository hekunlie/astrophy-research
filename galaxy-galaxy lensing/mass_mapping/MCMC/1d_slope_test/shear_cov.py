import numpy
from sys import path
path.append("E:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import matplotlib.pyplot as plt

seed = 1214
total_num = 10000
uni_num = 10000
rng = numpy.random.RandomState(seed)
fq = Fourier_Quad(6, 12)

bin_num = 8
bin_num2 = int(bin_num/2)
inverse = range(int(bin_num2 - 1), -1, -1)

# x = numpy.random.uniform(-8, 8, gal_num)
x_asy = rng.uniform(8, 8, int(total_num - uni_num))
x_uni = rng.uniform(-2, 8, uni_num)
x = numpy.zeros((total_num,))
x[:uni_num] = x_uni
x[uni_num:] = x_asy

a1, a2 = 0, -0.01
shear_slope_1d = a1 + a2*x

mg = rng.normal(shear_slope_1d, 0.3, total_num)

mnu = numpy.ones_like(mg)
mg_bins = fq.set_bin(mg, bin_num, 10)

img = Image_Plot()
img.subplots(1,1)
img.axs[0][0].scatter(x,mg,s=5)
img.show_img()

nx, ny = 41, 41
a1_range = numpy.linspace(-0.1, 0.1, nx)
a2_range = numpy.linspace(-0.1, 0.1, ny)

chisq = numpy.zeros((ny, nx))
chisq_new = numpy.zeros_like(chisq)
for i in range(nx*ny):
    m, n = divmod(i, nx)
    gh = a1_range[n] + a2_range[m]*x
    # y(a2)-axis, x(a1)-axis
    covs = numpy.cov(mg, mg - gh)
    chisq[m,n] = covs[1,1] - covs[0,0]
    # chisq_new[m,n] = fq.get_chisq_new(mg, 1, gh, mg_bins, bin_num2, inverse, 0, num_exp)

show_chisq_value = 20
norm = plt.Normalize(vmin=numpy.min(chisq), vmax=numpy.max(chisq))
norm_new = plt.Normalize(vmin=numpy.min(chisq_new), vmax=numpy.max(chisq_new))
norm_show = plt.Normalize(vmin=0, vmax=show_chisq_value)
cmap = plt.get_cmap('jet')

img = Image_Plot()
img.subplots(1,1)

for i in range(ny):
    for j in range(nx):
        # plot \chi (original) squared
        cl = cmap(norm(chisq[i, j]))
        img.axs[0][0].scatter(a1_range[j], a2_range[i], marker="s", color=cl, s=30)
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
# sm_new = plt.cm.ScalarMappable(cmap=cmap, norm=norm_new)
# sm_show = plt.cm.ScalarMappable(cmap=cmap, norm=norm_show)
sm._A = []
# sm_new._A = []
# sm_show._A = []
plt.colorbar(sm, ax=img.axs[0][0])

for i in range(2):
    img.set_label(0, 0, i, "$a_%d$" % (2 - i))


img.show_img()