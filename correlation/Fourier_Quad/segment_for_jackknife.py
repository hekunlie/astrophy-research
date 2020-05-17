import numpy
from sys import path
# path.append('%s/work/fourier_quad/'%my_home)
path.append("/home/hkli/work/mylib")
import matplotlib.pyplot as plt
import matplotlib
import h5py
from plot_tool import Image_Plot
from sklearn.cluster import KMeans


data_path = "/mnt/ddnfs/data_users/hkli/CFHT/correlation/cata/total.hdf5"


area_num = 4

h5f = h5py.File(data_path,"r")
for i in range(4):
    data = h5f["/w%d"%i][()]

    h5f.close()

ncen = 100

coord = data[:, 29:31]
km = kmeans_sample(coord, ncen, maxiter=100, tol=1.0e-5)
print(rank, km.converged, km.labels.size)



img= Image_Plot()
img.subplots(1,2)

sub_num = 300000
cmap = plt.cm.get_cmap('rainbow')

normalize = matplotlib.colors.Normalize(vmin=numpy.min(km.labels[:sub_num]), vmax=numpy.max(km.labels[:sub_num]))
# colors = [cmap(normalize(value)) for value in z]
colors = cmap(normalize(km.labels[:sub_num]))

img.axs[0][0].scatter(coord[:sub_num,0],coord[:sub_num,1],c=colors)
img.axs[0][1].hist2d(coord[:,0],coord[:,1],[200,200])

img.set_label(0,0,0,"Dec.")
img.set_label(0,0,1,"R.A.")

img.set_label(0,1,0,"Dec.")
img.set_label(0,1,1,"R.A.")

cax, _ = matplotlib.colorbar.make_axes(img.axs[0][0])
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)
img.save_img("w%d.png"%rank)
img.show_img()

comm.Barrier()

for i in range(numpros):
    if rank == i:
        if i == 0:
            h5f = h5py.File(data_path + "/labels.hdf5","w")
        else:
            h5f = h5py.File(data_path + "/labels.hdf5", "r")
        h5f["/w%d"%rank] = km.labels
        h5f.close()
    comm.Barrier()