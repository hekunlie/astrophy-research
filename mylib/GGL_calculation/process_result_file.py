from sys import path, argv
path.append("D:/Github/astrophy-research/mylib")
import shutil
import h5py
import numpy
import hk_tool_box
from hk_plot_tool import Image_Plot

data_path = "D:/TEMP/hjxu_cata/dr9/new"
file_name_rand = None #"result_rand_z_0.1_0.3_absz_-21.0_-22.0.hdf5"
file_name = "result_gal_z_0.1_0.3_absz_-20.0_-21.0.hdf5"


zhang_data = numpy.loadtxt(data_path + "/corr_2021.dat")
zhang_radius = zhang_data[:,5]
zhang_esd = zhang_data[:,2]
zhang_esd_err = zhang_data[:,3]

# for file_name in ["result_gal_z_0.1_0.3_absz_-20.0_-21.0.hdf5", "result_rand_z_0.1_0.3_absz_-20.0_-21.0.hdf5"]:
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

dilution_ratio = 1
if file_name_rand:
    h5f = h5py.File(data_path + "/%s" % file_name_rand, "r")
    rand_count_jkf = h5f["/count"][()]
    rand_radius_accum_jkf = h5f["/radius"][()]
    rand_theta_accum_jkf = h5f["/theta"][()]
    rand_len_count = h5f["/len_count"][()]
    h5f.close()
    rand_radius_jkf = rand_radius_accum_jkf/rand_count_jkf
    rand_theta_jkf = rand_theta_accum_jkf/rand_count_jkf

    dilution_ratio = (rand_count_jkf[jkf_num]/rand_len_count)/(count_jkf[jkf_num]/len_count)

    # print(rand_len_count, len_count)

delta_sigmat = delta_sigmat_jkf[:jkf_num].mean(axis=0)
delta_sigmat_err = delta_sigmat_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

delta_sigmax = delta_sigmax_jkf[:jkf_num].mean(axis=0)
delta_sigmax_err = delta_sigmax_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

radius = radius_jkf[:jkf_num].mean(axis=0)
radius_err = radius_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

theta = theta_jkf[:jkf_num].mean(axis=0)
theta_err = theta_jkf[:jkf_num].std(axis=0)*numpy.sqrt(jkf_num - 1)

pts_num = delta_sigmat.shape[0]
results = numpy.zeros((pts_num, 7))
results[:,0] = radius
results[:,1] = radius_err
results[:,2] = delta_sigmat
results[:,3] = delta_sigmat_err
results[:,4] = delta_sigmax
results[:,5] = delta_sigmax_err
results[:,6] = dilution_ratio

header = "Radius\tRadius_err\tDelta_Sigma\tDelta_Sigma_err\tDelta_Sigma_x\tDelta_Sigma_x_err\tDilution_ratio"
numpy.savetxt(data_path + "/%s.txt"%file_name, X=results, header=header)


img = Image_Plot(xpad=0.3,ypad=0.15)
img.subplots(1,2)
img.axs[0][0].errorbar(radius, delta_sigmat, delta_sigmat_err, fmt=" ", marker="s", capsize=3)
img.axs[0][0].errorbar(zhang_radius, zhang_esd, zhang_esd_err, fmt=" ", marker="s", capsize=3)
# img.axs[0][1].errorbar(radius, delta_sigmat*radius, delta_sigmat_err*radius, fmt=" ", marker="s", capsize=3)
# img.axs[0][2].errorbar(radius, delta_sigmax, delta_sigmax_err, fmt=" ", marker="s", capsize=3)
for i in range(1):
    img.axs[0][i].set_xscale("log")
    img.set_label(0,i,1, "Radius [Mpc/h]")
img.set_label(0,0,0,"$\Delta\Sigma$")
img.set_label(0,1,0,"$R\\times\Delta\Sigma$")
# img.set_label(0,2,0,"$\Delta\Sigma_{\\times}$")

img.axs[0][0].set_ylim(0.005,200)
img.axs[0][0].set_xlim(0.02, 100)
img.axs[0][0].set_yscale("log")
# img.save_img(data_path + "/%s.png"%file_name)
img.show_img()
