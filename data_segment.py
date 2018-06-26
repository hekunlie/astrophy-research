import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import time
from sys import argv
from subprocess import Popen
import shelve

fg_num = argv[1]

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_cut_path" in path:
        cut_path = path.split("=")[1]

data = numpy.load(result_path+"data_cache.npz")['arr_0']
cuts_num = 20
n_star = data[:, 3]
idx = n_star >= 20
area = data[:, 7]
a_idx = area >= 6

peak = data[:, 4][idx&a_idx]
flux = data[:, 5][idx&a_idx]
hflux = data[:, 6][idx&a_idx]
area = data[:, 7][idx&a_idx]
harea = data[:, 8][idx&a_idx]
flux2 = numpy.sqrt(data[:, 10][idx&a_idx])
flux_alt = numpy.sqrt(data[:, 11][idx&a_idx])

d_sort = numpy.sort(peak)
step = int(len(d_sort)/cuts_num)
peak_cut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(flux)
step = int(len(d_sort)/cuts_num)
flux_cut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(hflux)
step = int(len(d_sort)/cuts_num)
hflux_cut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(area)
step = int(len(d_sort)/cuts_num)
area_cut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(harea)
step = int(len(d_sort)/cuts_num)
harea_cut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(flux2)
step = int(len(d_sort)/cuts_num)
flux2_cut = [d_sort[i*step] for i in range(cuts_num)]

d_sort = numpy.sort(flux_alt)
step = int(len(d_sort)/cuts_num)
flux_alt_cut = [d_sort[i*step] for i in range(cuts_num)]

select = {"flux": (flux, flux_cut), "hflux": (hflux, hflux_cut), "area": (area, area_cut), "harea": (harea, harea_cut),
          "peak": (peak, peak_cut), "flux2":  (flux2, flux2_cut), "flux_alt": (flux_alt, flux_alt_cut)}
f = shelve.open(cut_path+"cut_dict")
f["dict"] = select
f.close()

# g1num = fg_num-10
# g2num = fg_num
# g1 = numpy.linspace(-0.005, 0.005, g1num)
# g2 = numpy.linspace(-0.0065, 0.0065, g2num)
# numpy.savez(result_path+"g_bin.npz", g1, g2)
# dg1 = g1[1]-g1[0]
# dg2 = g2[1]-g2[0]
#
# fg1 = data[:, 14]
# fg2 = data[:, 15]
#
# gs = [g1, g2]
# dgs = [dg1, dg2]
# fgs = [fg1, fg2]
# nums = [g1num, g2num]
# for i in range(2):
#     g = gs[i]
#     num = nums[i]
#     fg = fgs[i]
#     dg = dgs[i]
#     for j in range(num):
#         data_seg_name = result_path + "g%d_%d.npz"%(i+1,j)
#         idx1 = fg >= g[j] - dg/2.
#         idx2 = fg <= g[j] + dg/2.
#         data_seg = data[idx1&idx2]
#         numpy.savez(data_seg_name, data_seg)