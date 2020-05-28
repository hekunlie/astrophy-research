import h5py
import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot

h5f = h5py.File("test.hdf5","r")
img = h5f["/img"][()]
img_pow = h5f["/img_pow"][()]
img_pow_real = h5f["/img_pow_real"][()]
img_pow_imag = h5f["/img_pow_imag"][()]
img_pow_arg = h5f["/img_arg"][()]
h5f.close()

fq = Fourier_Quad(img.shape[0], 12)
img_pow_fq = fq.pow_spec(img)
img_pow_arg_fq = fq.pow_arg(img)
img_pow_real_fq, img_pow_imag_fq = fq.pow_real_imag(img)

datas = [[img, img_pow, img_pow_arg,img_pow_real,img_pow_imag],
         [img, img_pow_fq, img_pow_arg_fq,img_pow_real_fq,img_pow_imag_fq],
         [img - img, img_pow - img_pow_fq, img_pow_arg-img_pow_arg_fq,
          img_pow_real-img_pow_real_fq, img_pow_imag-img_pow_imag_fq]]

img = Image_Plot()
img.subplots(3,5)
for i in range(3):
    for j in range(5):
        fig = img.axs[i][j].imshow(datas[i][j])
        img.figure.colorbar(fig, ax=img.axs[i][j])
        img.del_ticks(i,j,[0,1])
img.show_img()