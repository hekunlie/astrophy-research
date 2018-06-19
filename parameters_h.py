import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
import numpy
import matplotlib.pyplot as plt
from Fourier_Quad import Fourier_Quad
import tool_box
import lsstetc
from mpi4py import MPI
import h5py


num = int(argv[1])*10000
size = int(argv[2])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

mag_s, mag_e = 20, 26.5
radius_s, radius_e = 0.7, 1.1

with open("%s/work/envs/envs1.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "select_parameter" in path:
        para_path = path.split("=")[1]

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
theta = rng.uniform(0, numpy.pi, num)
q = (1-e)/(1+e)
es = (1-q**2)/(1+q**2)
e1 = es*numpy.cos(2*theta)
e2 = es*numpy.sin(2*theta)
print("Rank: %3d, mean(e1): %10.6f, std: %.4f, mean(e2): %10.6f, std: %.4f, max: %.5f, %.5f"
      %(rank, numpy.mean(e1), numpy.std(e1), numpy.mean(e2), numpy.std(e2), numpy.max(e1), numpy.max(e2)))
f["/e1"] = e1
f["/e2"] = e2

# magnitude & flux

flux = numpy.zeros((num, 1))
arr = tool_box.mags_mock(num, mag_s, mag_e)
for i in range(num):
    flux[i] = prop.flux(arr[i])
f["/flux"] = flux[:,0]
f["/mag"] = arr

# galactic radius

rad = rng.uniform(radius_s, radius_e, num)
f["/radius"] = rad

# B/T ratio
bt = rng.normal(0, 0.1, 2*num)
bt.shape = (len(bt), 1)
idx1 = bt >= 0
idx2 = bt < 1
c = bt[idx1&idx2]
c.shape = (len(c), 1)
if len(c) > num:
    f_bt = c[:num]
elif len(c) < num:
    gap = num - len(c)
    plus = rng.normal(0, 0.1, 10*gap)
    plus = plus[plus>=0][:gap]
    plus.shape = (len(plus), 1)
    f_bt = numpy.row_stack((c, plus))
else:
    f_bt = c
f["/btr"] = f_bt
f.close()

pic = para_path + "/pic/e1e2ees_%d.png"%rank
plt.subplot(221)
plt.hist(e1, 100)
plt.subplot(222)
plt.hist(e2, 100)
plt.subplot(223)
plt.hist(e, 100)
plt.subplot(224)
plt.hist(es, 100)
plt.savefig(pic)
plt.close()
#
pic = para_path + "/pic/fmrb_%d.png"%rank
plt.subplot(221)
plt.hist(arr, 100)
plt.subplot(222)
plt.hist(flux[:,0], 100)
plt.subplot(223)
plt.hist(rad, 100)
plt.subplot(224)
plt.hist(f_bt, 100)
plt.savefig(pic)
plt.close()


