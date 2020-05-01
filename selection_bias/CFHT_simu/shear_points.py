import numpy
from sys import argv
import h5py
import matplotlib.pyplot as plt


# data_path, shear_num = argv[1], int(argv[2])

shear_num = 10
data_path = "D:/"
# g = numpy.linspace(0,0.05, shear_num)
# theta = numpy.random.uniform(0,2*numpy.pi,shear_num)
# cos_2theta = numpy.cos(2*theta)
# sin_2theta = numpy.sin(2*theta)
# g1 = g*cos_2theta
# g2 = g*sin_2theta
g1 = numpy.zeros((shear_num,))
g2 = numpy.zeros((shear_num,))
g1[:5] = numpy.linspace(-0.05, -0.005, 5)
g1[5:] = numpy.linspace(0.005, 0.05, 5)
g2[:5] = numpy.linspace(-0.05, -0.005, 5)
g2[5:] = numpy.linspace(0.005, 0.05, 5)
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