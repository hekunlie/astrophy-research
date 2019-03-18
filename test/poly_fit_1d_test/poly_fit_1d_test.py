import numpy
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
path.append('E:/Github/astrophy-research/my_lib/')
import tool_box
import h5py
from sys import argv
import matplotlib.pyplot as plt


num = int(argv[1])
order = int(argv[2])
percent = float(argv[3])


x = numpy.linspace(-5, 5, num)
fx = 0
coeff = numpy.random.uniform(-2,2,order+1)
for i in range(order+1):
    fx += coeff[i]*(x**i)
plt.plot(x, fx,label='origin')

print("Input coeff: ", coeff)
for i in range(num):
    fx[i] = fx[i] + numpy.random.normal(0, numpy.abs(percent*fx[i]),1)[0]

coeff,cov, fxy = tool_box.fit_1d(x, fx, order, "lsq")
coeff_s = tool_box.fit_1d(x, fx, order, "scipy")
# print(cov)
# print(fxy)
print("Lsq: ",coeff)
print("Scipy: ",coeff_s)
pic_title = "coeff: "
for i in range(order+1):
    pic_title += "%.4f, "%coeff[i]


fx_fit = 0
for i in range(order+1):
    fx_fit += coeff[i]*(x**i)
plt.scatter(x, fx, label="noise")
plt.plot(x, fx_fit,label="fit")
plt.title(pic_title)
plt.legend()
# plt.show()
plt.savefig("test.png")

f = h5py.File("test.hdf5", "w")
f["/x"] = x
f["/fx"] = fx
f.close()
