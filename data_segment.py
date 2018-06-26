import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import time
from sys import argv
from subprocess import Popen


fg_num = argv[1]

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_res_path" in path:
        result_path = path.split("=")[1]
data = numpy.load(result_path+"data_cache.npz")['arr_0']

g1num = fg_num-10
g2num = fg_num
g1 = numpy.linspace(-0.005, 0.005, g1num)
g2 = numpy.linspace(-0.0065, 0.0065, g2num)
numpy.savez(result_path+"g_bin.npz", g1,g2)
dg1 = g1[1]-g1[0]
dg2 = g2[1]-g2[0]

fg1 = data[:, 14]
fg2 = data[:, 15]

gs = [g1, g2]
dgs = [dg1, dg2]
fgs = [fg1, fg2]
nums = [g1num, g2num]
for i in range(2):
    g = gs[i]
    num = nums[i]
    fg = fgs[i]
    dg = dgs[i]
    for j in range(num):
        data_seg_name = result_path + "g%d_%d.npz"%(i+1,j)
        idx1 = fg >= g[j] - dg/2.
        idx2 = fg <= g[j] + dg/2.
        data_seg = data[idx1&idx2]
        numpy.savez(data_seg_name, data_seg)