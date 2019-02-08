import numpy
import os
from sys import path
path.append("E:/Github/astrophy-research/my_lib")
import tool_box
from astropy.io import fits
import matplotlib.pyplot as plt
import h5py

# gaps in the histogram!!!
edge = 34
pix_num = 200000
yd, xd = 4644-2*edge, 2112-2*edge
mask = numpy.zeros((yd, xd))

img = fits.open("E:/img.fits")[0].data
xyf = h5py.File("E:/ch_xy.hdf5")
xy = xyf["/data"].value
xyf.close()
h5f = h5py.File("E:/ch.hdf5")
data = h5f["/data"].value
h5f.close()

my, mx = numpy.mgrid[0:yd,0:xd]
my = my.flatten()
mx = mx.flatten()
ch_seqs = numpy.random.choice(numpy.arange(0,yd*xd), pix_num, replace=False)
img_f = img[edge:yd+edge,edge:xd+edge].flatten()[ch_seqs]
plt.hist(img_f,100)
plt.show()
fz_s = numpy.sort(img_f)
print(fz_s[5000])
bottom, upper = fz_s[int(pix_num*0.2)], fz_s[int(pix_num*0.8)]
idx_1 = fz_s >= bottom
idx_2 = fz_s < upper


ys = my[ch_seqs]
xs = mx[ch_seqs]

plt.hist(fz_s[idx_1&idx_2],100)
plt.show()


plt.imshow(img)
plt.show()
ch_list = []
for k in range(150000):
    i,j = divmod(xy[k,0], 2112-2*edge)
    mask[i,j] = 1
    ch_list.append(img[i,j])

ch_list = numpy.sort(numpy.array(ch_list))
print(ch_list.shape, ch_list[:10])
plt.hist(ch_list,100)
plt.show()
plt.imshow(mask)
plt.show()

print(data.shape)

for i in range(len(data)-1):
    if data[i,0] > data[i+1,0]:
        print(data[i,0],i)
a,b = 3321.15, 3370.5
idx1 = data<b
idx2 = data>=a

print(len(data[idx1&idx2]))
plt.hist(data[idx1&idx2], 50)
plt.show()
plt.imshow(data[50000:100000,0].reshape((200,250)))
plt.show()