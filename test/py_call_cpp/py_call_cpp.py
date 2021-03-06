import numpy
import ctypes
import numpy.ctypeslib as ctl
from numpy.ctypeslib import ndpointer
import time


lib = ctypes.cdll.LoadLibrary("libhist.so")


# test1 = lib.test1
# test1.restype = None
# test1.argtypes = [ctypes.c_double, ctypes.c_double]
# a,b = 10,15
# test1(a,b)

# this doesn't work, "segmentation fault"
# test2 = lib.test2
# test2.restype = None
# test2.argtypes = [ctypes.c_double, ctypes.c_double, ctypes.c_double]
# a,b = 10,15
# c = 1
# test2(a,b,c)
# print(c)



# test3 = lib.test3
# test3.restype = ctypes.c_double
# test3.argtypes = [ctypes.c_double, ctypes.c_double]
# a,b = 10,15
# c = test3(a, b)
# print(c)
#
#
# test4 = lib.test4
# test4.restype = ctypes.c_double
# test4.argtypes = [ctypes.c_int, ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous')]
# a = 1000
# b = numpy.ones((a, ), dtype=numpy.float64)
# c = test4(a, b)
# print(c)


# test = lib.test
# test.restype = None
# test.argtypes = [ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),ctypes.c_int,
#                  ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'), ctypes.c_int]
#
# num1, num2 = 1000, 10
# sum_val = ctypes.byref(ctypes.c_double())
# arr1 = numpy.ones((num1, ), dtype=numpy.float64)
# arr2 = numpy.zeros((num2, ), dtype=numpy.float64)
#
# test(arr1, num1, arr2, num2)
#
# print(arr2)



# locate = lib.locate
# locate.restype = ctypes.c_int
# locate.argtypes = [ctypes.c_double, ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'), ctypes.c_int]
#
# bin_num = 100
# a = numpy.random.uniform(-1000,1000,1)[0]
# bins = numpy.linspace(-1000,1000,bin_num+1)
#
# bin_tag = locate(a, bins, bin_num)
# print(bin_tag, bins[bin_tag],a, bins[bin_tag+1])



hist = lib.hist2d_fast
# hist_fast = lib.hist_fast

hist.restype = None
hist.argtypes = [ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctypes.c_int,
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctypes.c_int,
                 ctypes.c_int,
                 ctl.ndpointer(numpy.intc, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous')]

xbin_num = 500
ybin_num = 500
data_num = 3000000

x = numpy.random.uniform(-1000, 1000, data_num).astype(numpy.float64)
y = numpy.random.uniform(-1000, 1000, data_num).astype(numpy.float64)

xbins = numpy.linspace(-1000,1000,xbin_num+1, dtype=numpy.float64)
ybins = numpy.linspace(-1000,1000,ybin_num+1, dtype=numpy.float64)
count = numpy.zeros((ybin_num, xbin_num), dtype=numpy.intc)
grid_x = numpy.zeros_like(count, dtype=numpy.float64)
grid_y = numpy.zeros_like(count, dtype=numpy.float64)

t1 = time.time()

hist(x, y, data_num, xbins, ybins, xbin_num, ybin_num, count, grid_x, grid_y)
t2 = time.time()
num2d = numpy.histogram2d(y, x, [ybins, xbins])[0]
t3 = time.time()
diff = num2d - count
numpy.savez("./cache.npz", count, num2d, diff,grid_x,grid_y)
print(diff.min(), diff.max())
print(t2-t1, " sec")
print(t3-t2, " sec")





