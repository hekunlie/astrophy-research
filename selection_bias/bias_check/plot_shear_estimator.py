import matplotlib
# matplotlib.use("Agg")
from sys import path,argv
path.append("/home/hklee/work/mylib")
path.append("D:/GitHub/astrophy-research/mylib")
from plot_tool import Image_Plot
import numpy
import h5py
import tool_box
from Fourier_Quad import Fourier_Quad
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

data_path = argv[1]
shear_point = rank
data_type = argv[2]

fq = Fourier_Quad(12,124)

pic_path = data_path + "/diff_%d_%s.png"%(shear_point,data_type)
chi_pic_path = data_path + "/chi_%d_%s.png"%(shear_point,data_type)

h5f = h5py.File(data_path+"/shear.hdf5","r")
g1 = h5f["/g1"][()]
g2 = h5f["/g2"][()]
h5f.close()

print("g1: %.4f g2: %.4f"%(g1[shear_point],g2[shear_point]))

if data_type == "epsf":
    h5f = h5py.File(data_path+"/data_%d_noise_free_epsf.hdf5"%shear_point,"r")
    print("epsf")
else:
    h5f = h5py.File(data_path + "/data_%d_noise_free.hdf5" % shear_point, "r")
    print("cpsf")
data_noise_free_e = h5f["/data"][()]
h5f.close()

if data_type == "epsf":
    h5f = h5py.File(data_path+"/data_%d_noisy_cpp_epsf.hdf5"%shear_point,"r")
    print("epsf")
else:
    h5f = h5py.File(data_path + "/data_%d_noisy_cpp.hdf5" % shear_point, "r")
    print("cpsf")
data_noisy_e = h5f["/data"][()]
h5f.close()

scale = 10000

diff = (data_noise_free_e - data_noisy_e)

estimator = [diff[:,0],diff[:,2],diff[:,2]+diff[:,3],
             diff[:,1],diff[:,3],diff[:,2]-diff[:,3]]

est = ["G1","N","N+U",
       "G2","U","N-U"]

bin_num = 5000
img = Image_Plot()
img.subplots(2,3)

for i in range(6):
    m,n = divmod(i,3)

    img.axs[m][n].hist(estimator[i]/scale,bin_num)

    fit = tool_box.gaussnosie_fit(estimator[i]/scale,bin_num)
    img.axs[m][n].plot(fit[0],fit[1],c="C1",ls="--",linewidth=1.5)

    print("plot %s"%est[i])
    ys = img.axs[m][n].set_ylim()
    img.axs[m][n].plot([0,0],[ys[0],ys[1]],ls="--",alpha=0.5,c="k",linewidth=1)
    print(fit[4])
    text_str = "%s\n%.4f\n%.4f\n%.4f"%(est[i],fit[4][0],fit[4][1],fit[4][2])
    img.axs_text(m,n,0.8,0.1,text_str)

img.save_img(pic_path)

img = Image_Plot()
img.subplots(2,3)

for i in range(6):
    m,n = divmod(i,3)

    fq.find_shear(estimator[i],1,20,left=-3000, right=3000,fit_num=40,fig_ax=img.axs[m][n])

img.save_img(chi_pic_path)