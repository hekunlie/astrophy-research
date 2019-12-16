import numpy
from sys import path,argv
path.append("D:/Github/astrophy-research/mylib")
path.append("/home/hkli/work/mylib")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box
from astropy.io import fits
import time

num = int(argv[1])
size = int(argv[2])
cols = int(argv[3])
fq = Fourier_Quad(size,124)

img_row, img_col = divmod(num, cols)
img_stack_new = numpy.zeros((size * img_row, size * cols))


def stack_new(img, stamp, stamp_size, iy, ix):
    img[iy * stamp_size:(iy + 1) * stamp_size, ix * stamp_size:(ix + 1) * stamp_size] = stamp

for k in range(30):
    # t1 = time.time()
    # gals = []
    # for i in range(num):
    #     arr = numpy.zeros((size, size))+i+k
    #     gals.append(arr)
    # img_stack = fq.stack(gals, cols)
    # t2 = time.time()
    # print(t2-t1, img_stack.mean())
    # hdu = fits.PrimaryHDU(img_stack)
    # hdu.writeto("img.fits",overwrite=True)

    t3 = time.time()

    for i in range(num):
        iy, ix = divmod(i, cols)
        arr = numpy.zeros((size, size))+i+k
        fq.stack_new(img_stack_new, arr, iy, ix)

    t4 = time.time()
    print(t4-t3,img_stack_new.mean())
    # diff = img_stack_new - img_stack
    # print(diff.min(), diff.max())
    # hdu = fits.PrimaryHDU(img_stack_new)
    # hdu.writeto("img_new.fits",overwrite=True)