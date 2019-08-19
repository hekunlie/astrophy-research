import numpy
import matplotlib.pyplot as plt
from sys import path
path.append("E:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import h5py
from astropy.cosmology import FlatLambdaCDM
import tool_box

radius_bin = tool_box.set_bin_log(0.04,15, 13)

h = 0.7
C_0_hat = 2.99792458
H_0 = 100*h
coeff = 1000*C_0_hat/h

cosmos = FlatLambdaCDM(H_0, Om0=0.31)

fq = Fourier_Quad(12, 124)
filen = 8
result = numpy.zeros((5, filen))

rad_id = 6
area = 1
fore_zbin = 1
print(radius_bin[rad_id], radius_bin[rad_id+1])
data_path = "F:/works/GGL/%d/"%fore_zbin
cfht_data = numpy.loadtxt("E:/Github/astrophy-research/galaxy-galaxy lensing/GGL_calculation/lensing_low/data.dat")


h5f = h5py.File(data_path+"w%d/w_%d_%d.hdf5"%(area, area, fore_zbin), "r")
fore_redshift = h5f["/Z"].value
fore_dist = h5f["/DISTANCE"].value
print("Foreground Z: %.4f ~ %.4f. Mean dist: %.4f Mpc/h"
      %(fore_redshift.min(), fore_redshift.max(),fore_dist.mean()))

h5f = h5py.File(data_path+"w%d/cmass_result_w_%d.hdf5"%(area, area), "r")
result_f3 = h5f["/result"].value
x_f3 = h5f["/mean_dist"].value
h5f.close()


h5f = h5py.File(data_path + "w%d/radius_%d.hdf5" % (area,rad_id), "r")
data = h5f["/pair_data_0"].value
h5f.close()
print("The ESD: %.3f (%.3f)"%(result_f3[0,rad_id],result_f3[1,rad_id]))

crit_integ = data[:, 4] * 554.682135528

mgt = data[:, 0]
mgx = data[:, 1]
mnut = data[:, 2]
mnux = data[:, 3]
redshift = data[:,6]
dist_len = numpy.zeros_like(redshift)
for i in range(len(redshift)):
    dist_len[i] = cosmos.comoving_distance(redshift[i]).value
print("Source num: ", len(mgt))
print("Backgroud Z: %.4f ~ %.4f"%(redshift.min(), redshift.max()))


plt.hist(redshift, 20)
plt.show()
num_g = 6
result = numpy.zeros((5,num_g))

dz = 0.1
for i in range(num_g):
    z1 = redshift.min() + i*dz
    z2 = z1+dz
    idx_1 = redshift >= z1
    idx_2 = redshift < z2
    idx = idx_1 & idx_2

    # img = Image_Plot()
    # img.subplots(1, 2)
    left, right = -0.1, 0.1
    gt, gt_sig = fq.fmin_g_new(mgt[idx], mnut[idx], 8, left=left, right=right)[:2]
    gx, gx_sig = fq.fmin_g_new(mgx[idx], mnux[idx], 8, left=left, right=right)[:2]
    print("gamma_t: %.4f (%.4f).Num: %d.\nMead dist: %.4f Mpc/h. Z bin [%.4f, %.4f].\n"
          %(gt, gt_sig,idx.sum(),dist_len[idx].mean(),z1, z2))
    result[0, i] = gt
    result[1, i] = gt_sig
    result[2, i] = gx
    result[3, i] = gx_sig
    result[4, i] = dist_len[idx].mean()
    # img.show_img()
img = Image_Plot(fig_x=12, fig_y=8)
img.subplots(1, 1)
img.axs[0][0].errorbar(result[4],result[0], result[1], c="C4", label="$\gamma_t$")
img.set_label(0,0,0,"$\gamma_t$",size=15)
img.set_label(0,0,1,"$Dist_S\ [\\rm{Mpc \cdot h^{-1}}]$",size=15)
img.show_img()