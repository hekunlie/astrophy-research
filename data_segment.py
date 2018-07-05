import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import time
from sys import argv
from subprocess import Popen
import shelve
import shutil

fg_num = argv[1]

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_cut_path" in path:
        cut_path = path.split("=")[1]

data = numpy.load(result_path+"ori_data_cache.npz")['arr_0']
cuts_num = 20
n_star = data[:, 3]
idx = n_star >= 12
area = data[:, 7]
a_idx = area >= 6
s_data_path = result_path+"data_cache.npz"
s_data = data[idx&a_idx]
numpy.savez(s_data_path, data)

peak = s_data[:, 4]
flux = s_data[:, 5]
hflux = s_data[:, 6]
area = s_data[:, 7]
harea = s_data[:, 8]
flux2 = s_data[:, 10]
flux_alt = s_data[:, 11]

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

select = {"flux": flux_cut, "hflux": hflux_cut, "area": area_cut, "harea": harea_cut,
          "peak": peak_cut, "flux2":  flux2_cut, "flux_alt": flux_alt_cut}
f = shelve.open(cut_path+"cut_dict")
f["dict"] = select
f.close()

g1num = int(fg_num)-20
g2num = int(fg_num)

fg1_min, fg1_max = -0.005, 0.005
dg1 = (fg1_max - fg1_min)/g1num
g1 = numpy.array([fg1_min + (i-0.5)*dg1 for i in range(1,g1num+1)])


fg2_min, fg2_max = -0.0075, 0.0075
dg2 = (fg2_max - fg2_min)/g2num
g2 = numpy.array([fg2_min + (i-0.5)*dg2 for i in range(1,g2num+1)])
# g1 = numpy.linspace(-0.005, 0.005, g1num)
# g2 = numpy.linspace(-0.0075, 0.0075, g2num)
numpy.savez(result_path+"g_bin.npz", g1, g2)
# dg1 = g1[1]-g1[0]
# dg2 = g2[1]-g2[0]

fg1 = s_data[:, 14]
fg2 = s_data[:, 15]

gs = [g1, g2]
dgs = [dg1, dg2]
fgs = [fg1, fg2]
nums = [g1num, g2num]

seg_data_path = result_path + "data/"
shutil.rmtree(seg_data_path,ignore_errors=True)
os.mkdir(seg_data_path)

for i in range(2):
    g = gs[i]
    num = nums[i]
    fg = fgs[i]
    dg = dgs[i]
    for j in range(num):
        data_seg_name = result_path + "data/g%d_%d.npz"%(i+1,j)
        idx1 = fg >= g[j] - dg/2.
        idx2 = fg <= g[j] + dg/2.
        data_seg = s_data[idx1&idx2]
        numpy.savez(data_seg_name, data_seg)

# cut_fs = ["flux2", "flux_alt"]#, "peak", "flux", "area"]
# for cut_f in cut_fs:
#     cmd = "mpirun -np %s python sym_mc_plot_cfht.py %s"%(fg_num, cut_f)
#     a = Popen(cmd, shell=True)
#     a.wait()

# cmd = "mpirun -np %s python process_data.py 0 12 flux_alt 4"%fg_num
# a = Popen(cmd, shell=True)
# a.wait()