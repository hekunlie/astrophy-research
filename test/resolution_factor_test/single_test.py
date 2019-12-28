import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
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

total_path = argv[1]
shear_id = int(argv[2])
test_id = int(argv[3])
sex_filter = argv[4]

size = 50
num = 10000
fq = Fourier_Quad(size, 123)

chip_id,gal_label = divmod(test_id, num)
img_path = "%s/%d/gal_chip_%04d.fits"%(total_path, shear_id, chip_id)
img = fits.open(img_path)[0].data
gals = fq.segment(img)

rfactor_path = "%s/result/data/%s/rfactor_%d_new.hdf5"%(total_path, sex_filter, shear_id)
h5f_r = h5py.File(rfactor_path, "r")
rfactor_cpp = h5f_r["/data"][()]
h5f_r.close()

area_path = "%s/result/data/%s/area_%d.hdf5"%(total_path, sex_filter, shear_id)
h5f_a = h5py.File(area_path, "r")
area = h5f_a["/data"][()]
h5f_a.close()
print(img_path)
print(rfactor_path)
print(area_path)

print("Test on: shear %d, galaxy %d on chip %04d, %d gal, sex filter %s"%(shear_id, test_id, chip_id, gal_label, sex_filter))
print("Area: %d"%area[test_id])

gal_quad = tool_box.get_quad(gals[gal_label], size, numpy.sqrt(area[test_id]/numpy.pi))

print("CPP: %.6f, PY: %.6f"%(rfactor_cpp[test_id], gal_quad[0]))

img = Image_Plot(fig_x=4,fig_y=4)
img.subplots(1,2)
fig = img.axs[0][0].imshow(gals[gal_label])
img.figure.colorbar(fig, ax=img.axs[0][0])
fig = img.axs[0][1].imshow(gal_quad[2])
img.figure.colorbar(fig, ax=img.axs[0][1])
img.save_img("test.png")
