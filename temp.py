import matplotlib.pyplot as plt
import numpy
import time
import pandas
mag = numpy.loadtxt('E:/lsstmagsims')
mag = numpy.sort(mag)
print(numpy.min(mag),numpy.max(mag))
plt.hist(mag[0:10000],20)
plt.show()
plt.hist(mag[-10000:],20)
plt.show()