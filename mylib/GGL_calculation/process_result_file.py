from sys import path, argv
path.append("D:/Github/astrophy-research/mylib")
import shutil
import h5py
import numpy
import hk_tool_box
from hk_plot_tool import Image_Plot

data_path = "D:/TEMP"
file_name = "result_gal_z_0.1_0.3_absz_-20.0_-21.0.hdf5"
h5f = h5py.File(data_path + "/%s"%file_name,"r")
delta_sigmat_jkf = h5f["/delta_sigma_t"][()]
delta_sigmax_jkf = h5f["/delta_sigma_x"][()]
count_jkf = h5f["/count"][()]
radius_accum_jkf = h5f["/radius"][()]
theta_accum_jkf = h5f["/theta"][()]
len_count = h5f["/len_count"][()]
h5f.close()
radius_jkf = radius_accum_jkf/count_jkf
theta_jkf = theta_accum_jkf/count_jkf

jkf_num = count_jkf.shape[0] - 1
print("%d Jackknife"%jkf_num)

delta_sigmat = delta_sigmat_jkf[jkf_num]
delta_sigmat_err = delta_sigmat_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

delta_sigmax = delta_sigmax_jkf[jkf_num]
delta_sigmax_err = delta_sigmax_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

radius = radius_jkf[jkf_num]
radius_err = radius_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

theta = theta_jkf[jkf_num]
theta_err = theta_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

pts_num = delta_sigmat.shape[0]
results = numpy.zeros((pts_num, 6))
results[:,0] = radius
results[:,1] = radius_err
results[:,2] = delta_sigmat
results[:,3] = delta_sigmat_err
results[:,4] = delta_sigmax
results[:,5] = delta_sigmax_err

header = "Radius\tRadius_err\tDelta_Sigma\tDelta_Sigma_err\tDelta_Sigma_x\tDelta_Sigma_x_err"
numpy.savetxt("D:/TEMP/result_gal_z_0.1_0.3_absz_-20.0_-21.0.txt", X=results, header=header)

img = Image_Plot(xpad=0.3,ypad=0.15)
img.subplots(1,2)
img.axs[0][0].errorbar(radius, delta_sigmat, delta_sigmat_err, fmt=" ", marker="s", capsize=3)
img.axs[0][1].errorbar(radius, delta_sigmax, delta_sigmax_err, fmt=" ", marker="s", capsize=3)
for i in range(2):
    img.axs[0][i].set_xscale("log")
    img.set_label(0,i,1, "Radius [Mpc/h]")
img.axs[0][0].set_yscale("log")
img.save_img("D:/TEMP/result.png")
img.show_img()
