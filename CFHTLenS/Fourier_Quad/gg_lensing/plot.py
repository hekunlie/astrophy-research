from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
import h5py
from plot_tool import Image_Plot
import numpy

img = Image_Plot(xpad=0.25,ypad=0.2)
img.subplots(1,2)
nstar = [12,16,20,24]
names = ["result_w1", "result_w3"]
for i in range(2):
    for j in range(4):
        data_path = "E:/works/Galaxy-Galaxy_lensing/%s_%d.hdf5"%(names[i], j)
        h5f = h5py.File(data_path,"r")
        theta = h5f["/radius"][()][200]
        count = h5f["/count"][()][200]
        delta_t = h5f["/delta_t"][()][200]
        delta_t_err = h5f["/delta_t_err"][()][200]
        delta_x = h5f["/delta_x"][()][200]
        delta_x_err = h5f["/delta_x_err"][()][200]
        h5f.close()


        img.axs[0][i].errorbar(theta/count, delta_t, delta_t_err,capsize=3,marker="s", fmt=" ", label="Nstar >= %d"%nstar[j])
        # img.axs[1][i].errorbar(theta/count, delta_x, delta_x_err, capsize=3,marker="s",label="Nstar >= %d"%nstar[j])
        # print(i)
for i in range(2):
    img.set_label(0,i,0,"$\Delta\Sigma(R)$ [$h\cdot M_\odot\cdot pc^{-2}$]")
    img.set_label(0,i,1,"$R$ [Mpc/h]")
    img.axs[0][i].legend()
    img.axs[0][i].set_xscale("log")
    img.axs[0][i].set_yscale("log")
    img.axs[0][i].set_ylim(0.1,180)
    # for j in range(2):
    #     img.axs[i][j].legend()
    #     img.axs[i][j].set_xscale("log")
    #     img.axs[0][j].set_yscale("log")
img.save_img("D:/TEMP/ggl.png")
img.show_img()