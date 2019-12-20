import numpy
from sys import argv
import h5py
import matplotlib.pyplot as plt


# data_path, shear_num = argv[1], int(argv[2])

shear_num = 10
data_path = "D:/"

g1 = numpy.linspace(-0.04, 0.04, shear_num)
g2 = numpy.linspace(-0.04, 0.04, shear_num)
print(g1)
numpy.random.shuffle(g1)
print(g1)
numpy.random.shuffle(g2)
g = numpy.zeros((2*shear_num,))
g[:shear_num] = g1
g[shear_num:] = g2
plt.subplot(131)
plt.scatter(g1,g2)
plt.subplot(132)
plt.scatter(g1,g1)
plt.subplot(133)
plt.scatter(g2,g2)
plt.show()
h5f = h5py.File(data_path+"shear.hdf5","w")
h5f["/g1"] = g1
h5f["/g2"] = g2
h5f.close()

numpy.savetxt(data_path+"shear.dat",g)