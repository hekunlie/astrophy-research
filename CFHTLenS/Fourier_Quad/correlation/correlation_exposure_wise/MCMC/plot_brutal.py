from sys import path, argv
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import h5py
import numpy



data_path = "E:/works/correlation/CFHT/cut_2.5/deblend_1"

h5f = h5py.File(data_path + "/chisq_diff_expo.hdf5","r")
chisq = h5f["/chisq"][()]
omg_b = h5f["/omega_bm_bin"][()]
omg_c = h5f["/omega_cm_bin"][()]
sig8 = h5f["/sigma8"][()]

x, y = numpy.zeros_like(chisq[0]), numpy.zeros_like(chisq[0])
for i in range(len(sig8)):
    x[:,i] = sig8[i]
for i in range(len(omg_c)):
    y[i] = omg_c[i]

img = Image_Plot(fig_x=5,fig_y=4,xpad=0.25,ypad=0.25)
img.subplots(2,5)
for i in range(2):
    for j in range(5):
        chisq_i = -chisq[int(i*5 + j)]

        idx = chisq_i < 70

        xi = x[idx].flatten()
        yi = y[idx].flatten()
        zi = chisq_i[idx].flatten()

        idx = zi == zi.min()
        xi_min, yi_min = xi[idx], yi[idx]
        print(xi_min, yi_min)
        img.scatter_pts(i,j, xi,yi,numpy.log10(zi),color_map='jet',pts_size=5)
        img.axs[i][j].scatter(xi_min, yi_min, s=15,c="r", marker="p",label="$\chi^2_{min}=%.3f$"%zi.min())
        img.axs[i][j].legend()
        # fig = img.axs[i][j].imshow(chisq_i,cmap="jet")
        # img.figure.colorbar(fig,ax=img.axs[i][j],)
        img.set_label(i,j,0,"$\Omega_m$")
        img.set_label(i,j,1,"$\sigma_8$")
        img.axs[i][j].set_title("$\Omega_b=%.2f$"%omg_b[int(i*5 + j)])
# img.save_img(data_path +"/brutal.pdf")
img.show_img()
