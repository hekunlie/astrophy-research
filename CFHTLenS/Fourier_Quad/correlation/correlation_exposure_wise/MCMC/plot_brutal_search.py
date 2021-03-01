import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
import h5py
from plot_tool import Image_Plot



chisa_name = "chisq_diff_expo"
data_path = "E:/works/correlation/CFHT/cut_2.5/smooth/zbin_111111"

h5f = h5py.File("%s/%s.hdf5"%(data_path,chisa_name), "r")
chisqs = h5f["/chisq"][()]
omega_bm_bin = h5f["/omega_bm_bin"][()]
omega_cm_bin = h5f["/omega_cm_bin"][()]
sigma8_bin = h5f["/sigma8"][()]
h5f.close()
print(chisqs.shape)
print(chisqs.min(), chisqs.max())
omega_cm_num = omega_cm_bin.shape[0]
omega_bm_num = omega_bm_bin.shape[0]
sigma8_num = sigma8_bin.shape[0]

my,mx = numpy.zeros((omega_cm_num, sigma8_num)),numpy.zeros((omega_cm_num, sigma8_num))
for i in range(omega_cm_num):
    my[i] = omega_cm_bin[i]
for j in range(sigma8_num):
    mx[:,j] = sigma8_bin[j]


img = Image_Plot(fig_x=4,fig_y=3,xpad=0.25,ypad=0.28)
img.subplots(2,5)
for i in range(2):
    for j in range(5):
        tag = i*2 + j

        idx = chisqs[tag] > -200
        print(idx.sum())
        x = mx[idx].flatten()
        y = my[idx].flatten()
        z = chisqs[tag][idx].flatten()

        img.scatter_pts(i,j,x,y,z,color_map="jet",pts_size=18)

        img.axs[i][j].set_title("$\Omega_b=%.3f$"%omega_bm_bin[tag])
        img.set_label(i,j,0,"$\Omega_m$")
        img.set_label(i,j,1,"$\sigma_8$")
img.save_img("%s/%s.png"%(data_path, chisa_name))
img.show_img()