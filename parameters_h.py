from sys import path
path.append("/home/hkli/work/fourier_quad")
import numpy
import matplotlib.pyplot as plt
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
import time
import tool_box
import lsstetc
from mpi4py import MPI
import h5py


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "parameter" in path:
        para_path = path.split("=")[1]


num = 10000000
size = 84

seed = rank*4321554 + int(numpy.random.randint(1, 1256542344, 1)[0])
rng = numpy.random.RandomState(seed)
fq = Fourier_Quad(size, seed)
prop = lsstetc.ETC(band='r', pixel_scale=0.2, stamp_size=size, nvisits=180)
path = para_path+'para_%d.hdf5'%rank
f = h5py.File(path,"w")

# def ran(start, end, num):
#     return numpy.random.uniform(start, end, num)
#
#
# g1 = numpy.append(numpy.append(ran(-0.02, 0, 4), ran(0, 0.021, 3)),numpy.append(ran(-0.02, 0, 3), ran(0, 0.021, 4)))
# g2 = numpy.append(numpy.append(ran(-0.02, 0, 3), ran(0, 0.021, 4)),numpy.append(ran(0, 0.021, 3), ran(-0.02, 0, 4)))
# # numpy.random.shuffle(g1)
# # numpy.random.shuffle(g2)
# plt.subplot(131)
# plt.scatter(g1,g2)
# plt.subplot(132)
# plt.scatter(g1,g1)
# plt.subplot(133)
# plt.scatter(g2,g2)
# plt.show()
# numpy.savez('E:/selection_bias/parameters/shear.npz', g1, g2)
# numpy.savetxt('E:/selection_bias/parameters/shear.dat',numpy.append(g1,g2))


# ellipticity

e = tool_box.ellip_mock(num, seed)
theta = numpy.random.uniform(0, numpy.pi, num)
q = (1-e)/(1+e)
es = (1-q**2)/(1+q**2)
e1 = es*numpy.cos(2*theta)
e2 = es*numpy.sin(2*theta)
print("Rank: %3d, mean(e1): %10.6f, std: %.4f, mean(e2): %10.6f, std: %.4f, max: %.5f, %.5f"
      %(rank, numpy.mean(e1), numpy.std(e1), numpy.mean(e2), numpy.std(e2), numpy.max(e1), numpy.max(e2)))
f["/e1"] = e1
f["/e2"] = e2

# # magnitude & flux

flux = numpy.zeros((num, 1))
arr = tool_box.mags_mock(num, 18.5, 25)
for i in range(num):
    flux[i] = prop.flux(arr[i])
f["/flux"] = flux[:,0]



# galactic radius

rad = numpy.random.uniform(0.8, 1.8, num)
f["/radius"] = rad

f.close()
