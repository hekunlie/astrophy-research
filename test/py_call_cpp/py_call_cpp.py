import numpy
import ctypes
import numpy.ctypeslib as ctl
from numpy.ctypeslib import ndpointer


lib = ctypes.cdll.LoadLibrary("libhist.so")

hist = lib.hist2d
# hist_fast = lib.hist_fast

hist.restype = None
hist.argtypes = [ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctypes.c_int,
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctypes.c_int,
                 ctypes.c_int,
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous'),
                 ctl.ndpointer(numpy.float64, flags='aligned, c_contiguous')]

xbin_num = 1000
ybin_num = 800
data_num = 30000

x = numpy.random.uniform(-1000, 1000, data_num).astype(numpy.float64)
y = numpy.random.uniform(-1000, 1000, data_num).astype(numpy.float64)

xbins = numpy.linspace(-1000,1000,xbin_num+1, dtype=numpy.float64)
ybins = numpy.linspace(-1000,1000,ybin_num+1, dtype=numpy.float64)
count = numpy.zeros((ybin_num, xbin_num), dtype=numpy.float64)
grid_x = numpy.zeros_like(count, dtype=numpy.float64)
grid_y = numpy.zeros_like(count, dtype=numpy.float64)

hist(x, y, data_num, xbins, ybins, xbin_num, ybin_num, count, grid_x, grid_y)

num2d = numpy.histogram2d(y, x, [ybins, xbins])[0]
diff = num2d - count
print(diff.min(), diff.max())




