from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from plot_tool import Image_Plot
import numpy
import tool_box


data = numpy.load("D:/TEMP/chain_diff_expo_autocorr_thin_20000_steps.npz")["arr_0"]

As = data[:,0]/10**9
omega_cm0 = data[:,1]*(1-data[:,2])
omega_bm0 = data[:,1]*data[:,2]

data_s8 = numpy.load("D:/TEMP/chain_diff_expo_autocorr_thin_20000_steps_s8.npz")["arr_0"]

sig8 = data_s8[:,0]

diff_cm = omega_cm0 - data_s8[:,1]
diff_bm = omega_bm0 - data_s8[:,2]

print(diff_cm.min(), diff_cm.max())
print(diff_bm.min(), diff_bm.max())

# img = Image_Plot(fig_x=5,fig_y=4, xpad=0.2,ypad=0.25)
# img.subplots(2,2)
#
# img.axs[0][0].hist(sig8,100)
# img.set_label(0,0,1,"$\sigma_8$")
#
# img.axs[0][1].hist(As,100)
# img.set_label(0,1,1,"$A_s$")
#
# img.axs[1][0].hist(omega_cm0,100)
# img.set_label(1,0,1,"$\Omega_m$")
#
# img.axs[1][1].hist(omega_bm0,100)
# img.set_label(1,1,1,"$\Omega_b$")
# img.save_img("D:/TEMP/para_pdf.png")
# img.show_img()

# lower, upper = 0.538,0.54
# xs = (lower - (upper-lower)*0.2,upper + (upper-lower)*0.2)
# idx1 = sig8 <= upper
# idx2 = sig8 >= lower
# idx = idx1&idx2
# print(idx.sum())
#
# bin_num = 30
# img = Image_Plot(fig_x=5,fig_y=4, xpad=0.25,ypad=0.25)
# img.subplots(1,3)
#
# img.axs[0][0].scatter(sig8[idx],As[idx]*10**9)
# img.set_label(0,0,1,"$\sigma_8$")
# img.set_label(0,0,0,"$A_s$")
# img.axs[0][0].set_xlim(xs)
#
# img.axs[0][1].scatter(sig8[idx],omega_cm0[idx])
# img.set_label(0,1,1,"$\sigma_8$")
# img.set_label(0,1,0,"$\Omega_m$")
# img.axs[0][1].set_xlim(xs)
#
# img.axs[0][2].scatter(sig8[idx],omega_bm0[idx])
# img.set_label(0,2,1,"$\sigma_8$")
# img.set_label(0,2,0,"$\Omega_b$")
# img.axs[0][2].set_xlim(xs)
#
# img.save_img("D:/TEMP/para_relation.png")
# img.show_img()

sig8_sort = numpy.sort(sig8)
num = len(sig8)
mid = int(num*0.5)-1
sig8_best = numpy.percentile(sig8,50)

diff = si


