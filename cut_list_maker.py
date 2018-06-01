import matplotlib
matplotlib.use("Agg")
import numpy
import os
from sys import path, argv
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import shelve


with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_res_path" in path:
        result_path = path.split("=")[1]
    elif "cfht_pic_path" in path:
        pic_path = path.split("=")[1]


data_cache = result_path + "data_cache.npz"
data = numpy.load(data_cache)['arr_0']
cuts_num = int(float(argv[1]))
n_star = data[:, 3]
idx = n_star >= int(float(argv[2]))

nsig = data[:, 4][idx]
flux = data[:, 5][idx]
hflux = data[:, 6][idx]
area = data[:, 7][idx]
harea = data[:, 8][idx]
fsnr = numpy.sqrt(data[:, 10][idx])
fsnr_f = numpy.sqrt(data[:, 11][idx])
snr08 = numpy.sqrt(flux)/nsig

d_sort = numpy.sort(flux)
step = int(len(d_sort)/cuts_num)
fcut = [d_sort[i*step] for i in range(cuts_num)]
dict_path = result_path + "flux.npz"
numpy.savez(dict_path, fcut)

d_sort = numpy.sort(snr08)
step = int(len(d_sort)/cuts_num)
scut = [d_sort[i*step] for i in range(cuts_num)]
dict_path = result_path + "snr.npz"
numpy.savez(dict_path, scut)

d_sort = numpy.sort(fsnr)
step = int(len(d_sort)/cuts_num)
fsnrcut = [d_sort[i*step] for i in range(cuts_num)]
dict_path = result_path + "fsnr.npz"
numpy.savez(dict_path, fsnrcut)

d_sort = numpy.sort(fsnr_f)
step = int(len(d_sort)/cuts_num)
ffsnrcut = [d_sort[i*step] for i in range(cuts_num)]
dict_path = result_path + "fsnr_f.npz"
numpy.savez(dict_path, ffsnrcut)


