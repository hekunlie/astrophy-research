import numpy
import h5py
from subprocess import Popen
from sys import argv
import matplotlib.pyplot as plt

def get_chi_2d(num_arr, nx, ny):
    # arr_1 = num_arr[0:ny, 0:nx][:,range(ny-1,-1,-1)]
    # arr_2 = num_arr[0:ny, nx:2*nx]
    # arr_3 = num_arr[ny:2*ny, 0:nx][range(ny-1,-1,-1)][:,range(nx-1,-1,-1)]
    # arr_4 = num_arr[ny:2*ny, nx:2*nx][range(ny-1,-1,-1)]
    # chi_sq = 0.5 * numpy.sum(((arr_2 + arr_3 - arr_1 - arr_4) ** 2) / (arr_1 + arr_2 + arr_3 + arr_4))

    arr_1 = num_arr[ny:,nx:]
    arr_2 = num_arr[ny:, :nx][:,range(ny-1,-1,-1)]
    arr_3 = num_arr[:ny, :nx][range(ny-1,-1,-1)][:,range(ny-1,-1,-1)]
    arr_4 = num_arr[:ny, nx:][range(ny-1,-1,-1)]
    chi_sq = 0.5 * numpy.sum(((arr_1 + arr_3 - arr_2 - arr_4) ** 2) / (arr_1 + arr_2 + arr_3 + arr_4))

    # print(numpy.sum((arr_2 + arr_3 - arr_1 - arr_4) ** 2), numpy.sum(arr_1 + arr_2 + arr_3 + arr_4))
    return chi_sq

nx = int(argv[1])
num = numpy.random.normal(0,5, (nx,nx))
# my, mx = numpy.mgrid[0:nx,0:nx]
# num = numpy.exp(-((my-nx/2+0.5)**2 + (mx-nx/2+0.5)**2)/2/(0.2*nx)**2)*1000
# plt.imshow(num)
# plt.show()
h5f = h5py.File("test.hdf5","w")
h5f["/data"] = num
h5f.close()
# print(num)
chisq = get_chi_2d(num, int(nx/2), int(nx/2))
print(chisq)

cmd = "mpirun -np 1 ./get_corr %d"%nx
a = Popen(cmd, shell="True")
a.wait()


