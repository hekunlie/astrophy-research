import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from plot_tool import Image_Plot
import component_fit
import matplotlib.pyplot as plt

# data
pts_num = 5000000
g1, g2 = 0.1, -0.1
rng = numpy.random.RandomState(333)
mg = rng.multivariate_normal([g1,g2],[[0.8,0],[0,0.8]],pts_num)

# noise free data
mg1 = mg[:,0]#rng.normal(g1, 0.8, pts_num)
mg1_sym = mg1 - g1
mg2 = mg[:,1]#rng.normal(g2, 0.8, pts_num)
mg2_sym = mg2 - g2

# the scatter
para = [0.2, 0, 0.65]
sigma_g1 = para[0] + para[1]*mg1 + para[2]*mg1**2
sigma_g2 = para[0] + para[1]*mg2 + para[2]*mg2**2
# scattered data
mg1_scatter = numpy.random.normal(0, sigma_g1, pts_num)
mg1_ny = mg1_scatter + mg1
mg1_ny_sym = mg1_ny - g1
mg2_scatter = numpy.random.normal(0, sigma_g2, pts_num)
mg2_ny = mg2_scatter + mg2
mg2_ny_sym = mg2_ny - g2

# show the scatter
bin_num = 40
bound_nf = numpy.abs(mg1).max()/2
bins_nf = numpy.linspace(-bound_nf, bound_nf, bin_num+1)
bins_mid = (bins_nf[1:] + bins_nf[:-1])/2

fx = para[0] + para[1]*bins_mid + para[2]*bins_mid**2

img = Image_Plot(xpad=0.25,ypad=0.2)
img.subplots(1,2)
num = img.axs[0][0].hist(mg1, bins_nf)[0]
img.axs[0][0].errorbar(bins_mid, num, xerr=fx,capsize=3)

num = img.axs[0][1].hist(mg1-g1, bins_nf-g1)[0]
img.axs[0][1].errorbar(bins_mid-g1, num, xerr=fx,capsize=3)

for i in range(2):
    ys = img.axs[0][i].set_ylim()
    img.axs[0][i].plot([0,0], [ys[0],ys[1]], c="gray", ls="--")
    img.set_sci_ax(0,i,0)
    img.set_label(0,i,0,"Num")
    img.set_label(0,i,1,"G")
img.show_img()


# bin difference before and after PDF
bin_num = 200
num_ave_after_pdf = numpy.zeros((2,bin_num))

img = Image_Plot(xpad=0.25,ypad=0.2)
img.subplots(2,2)

# noise free data
bound = numpy.abs(mg1).max()/1.5
bins = numpy.linspace(-bound, bound, bin_num+1)
bins_mid = (bins[1:] + bins[:-1])/2

num = numpy.histogram(mg1, bins)[0]
num_sym = numpy.histogram(mg1_sym, bins)[0]

for i in range(int(bin_num/2)):
    num_ave_after_pdf[0,i] = (num_sym[i] + num_sym[bin_num-1-i])/2
    num_ave_after_pdf[0,bin_num-1-i] = (num_sym[i] + num_sym[bin_num-1-i])/2


img.axs[0][0].plot(bins_mid, num - num_sym, label="diff before & after PDF")
img.axs[0][0].plot(bins_mid, num, label="G1")
img.axs[0][0].plot(bins_mid, num_sym, label="G1-PDF")
img.axs[0][0].plot(bins_mid, num_sym - num_ave_after_pdf[0], label="G1-dipole after PDF")


# noisy data
# bound = numpy.abs(mg_ny).max()/10
# bins = numpy.linspace(-bound, bound, bin_num+1)
# bins_mid = (bins[1:] + bins[:-1])/2

num_ny = numpy.histogram(mg1_ny, bins)[0]
num_ny_sym = numpy.histogram(mg1_ny_sym, bins)[0]

for i in range(int(bin_num/2)):
    num_ave_after_pdf[1,i] = (num_ny_sym[i] + num_ny_sym[bin_num-1-i])/2
    num_ave_after_pdf[1,bin_num-1-i] = (num_ny_sym[i] + num_ny_sym[bin_num-1-i])/2

img.axs[0][1].plot(bins_mid, num_ny - num_ny_sym, label="diff before & after PDF")
img.axs[0][1].plot(bins_mid, num_ny, label="G1")
img.axs[0][1].plot(bins_mid, num_ny_sym, label="G1-PDF")
img.axs[0][1].plot(bins_mid, num_ny_sym - num_ave_after_pdf[1], label="G1-dipole after PDF")
for i in range(2):
    for j in range(2):
        ys = img.axs[j][i].set_ylim()
        xs = img.axs[j][i].set_xlim()
        img.axs[j][i].plot([0,0], [ys[0],ys[1]], c="gray", ls="--")
        img.axs[j][i].plot([xs[0],xs[1]], [0,0], c="gray", ls="--")
        img.set_label(j,i,0,"Num")
        img.set_label(j,i,1,"G")
        img.set_sci_ax(j,i,0)
        img.axs[j][i].legend()
img.show_img()

# 2D
img = Image_Plot(xpad=0.25, ypad=0.25)
img.subplots(2,2)

# noise free data
num, xgrid, ygrid = component_fit.get_2dhist(mg1, mg2, bins)[:3]

num = num.flatten()
xgrid = xgrid.flatten()
ygrid = ygrid.flatten()

norm = plt.Normalize(vmin=numpy.min(num), vmax=numpy.max(num))
cmap = plt.get_cmap('jet')
cl = cmap(norm(num))

img.axs[0][0].scatter(xgrid,ygrid,color=cl,s=2,marker="s")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
img.figure.colorbar(sm, ax=img.axs[0][0])

# noise free after PDF
num, xgrid, ygrid, radius = component_fit.get_2dhist(mg1_sym, mg2_sym, bins)[:4]
dipole = component_fit.get_dipole(num,radius,int(bin_num/3))[0]

dipole = dipole.flatten()
xgrid = xgrid.flatten()
ygrid = ygrid.flatten()

norm = plt.Normalize(vmin=numpy.min(dipole), vmax=numpy.max(dipole))
cmap = plt.get_cmap('jet')
cl = cmap(norm(dipole))

img.axs[1][0].scatter(xgrid,ygrid,color=cl,s=2,marker="s")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
img.figure.colorbar(sm, ax=img.axs[1][0])


# noisy data
num, xgrid, ygrid, radius = component_fit.get_2dhist(mg1_ny, mg2_ny, bins)[:4]
dipole = component_fit.get_dipole(num,radius,int(bin_num/3))[0]

num = num.flatten()
xgrid = xgrid.flatten()
ygrid = ygrid.flatten()

norm = plt.Normalize(vmin=numpy.min(num), vmax=numpy.max(num))
cmap = plt.get_cmap('jet')
cl = cmap(norm(num))

img.axs[0][1].scatter(xgrid,ygrid,color=cl,s=2,marker="s")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
img.figure.colorbar(sm, ax=img.axs[0][1])

# noisy data after PDF
num, xgrid, ygrid, radius = component_fit.get_2dhist(mg1_ny_sym, mg2_ny_sym, bins)[:4]
dipole = component_fit.get_dipole(num,radius,int(bin_num/3))[0]

dipole = dipole.flatten()
xgrid = xgrid.flatten()
ygrid = ygrid.flatten()

norm = plt.Normalize(vmin=numpy.min(dipole), vmax=numpy.max(dipole))
cmap = plt.get_cmap('jet')
cl = cmap(norm(dipole))

img.axs[1][1].scatter(xgrid,ygrid,color=cl,s=2,marker="s")
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm._A = []
img.figure.colorbar(sm, ax=img.axs[1][1])

for i in range(2):
    for j in range(2):
        ys = img.axs[i][j].set_ylim()
        xs = img.axs[i][j].set_xlim()
        img.axs[i][j].plot([xs[0],xs[1]],[0,0],c="gray",alpha=0.4,ls="--")
        img.axs[i][j].plot([0,0],[ys[0],ys[1]],c="gray",alpha=0.4,ls="--")
        img.set_label(i,j,0,"G1")
        img.set_label(i,j,1,"G2")
img.show_img()