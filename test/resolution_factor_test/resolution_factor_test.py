import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from astropy.io import fits
import time
from plot_tool import Image_Plot


t1 = time.time()




size = 64
num = 10000
fq = Fourier_Quad(size, 123)

img = fits.open("/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/pts_dimmer/1/gal_chip_0000.fits")[0].data
h5f_sex = h5py.File("/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/pts_dimmer/result/data/sex2_1.5/sex_1.hdf5","r")
sex_data = h5f_sex["/data"].value
h5f_sex.close()

gals = fq.segment(img)

quads = numpy.zeros((num,))

for i in range(len(gals)):
    eff_radius = sex_data[i,3]/numpy.pi
    if eff_radius > 0:
        quads[i] = tool_box.get_quad(gals[i], eff_radius, size)
    # if i < 20:
    #     print(sex_data[i,3],eff_radius, quads[i])

h5f = h5py.File("quad_size.hdf5", "r")
quad_cpp = h5f["/data"].value[:,0]
h5f.close()

# img = Image_Plot()
# img.subplots(1,1)
# img.axs[0][0].hist(quads,100)
# img.save_img("quad_hist.png")

diff = quads - quad_cpp
t2 =time.time()

print(diff.max(), diff.min(), t2-t1)
