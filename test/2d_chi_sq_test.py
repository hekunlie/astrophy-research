import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
import numpy
import h5py

def chi2d(num_arr):
    y,x = num_arr.shape
    ny, nx = int(y*0.5), int(x*0.5)
    arr_1 = num_arr[0:ny, 0:nx][:, range(ny - 1, -1, -1)]
    arr_2 = num_arr[0:ny, nx:2 * nx]
    arr_3 = num_arr[ny:2 * ny, 0:nx][range(ny - 1, -1, -1)][:, range(nx - 1, -1, -1)]
    arr_4 = num_arr[ny:2 * ny, nx:2 * nx][range(ny - 1, -1, -1)]
    return 0.5 * numpy.sum(((arr_2 + arr_3 - arr_1 - arr_4) ** 2) / (arr_1 + arr_2 + arr_3 + arr_4))

size = int(argv[1])

arr = numpy.random.randint(1,1000000, size*size,dtype="uint64").reshape((size, size))
h5f = h5py.File("hist_arr.hdf5","w")
h5f["/data"] = arr
h5f.close()

chi_sq = chi2d(arr)
print(chi_sq)

for i in range(size):
    s = ""
    for j in range(size):
        s += "%d, "%arr[i,j]
    print(s)
exit()