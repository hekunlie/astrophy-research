import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import time
import tool_box
from Fourier_Quad import *
from sys import argv
import numpy
import lsstetc
import h5py


ch, snr_s, scale, del_bin = argv[1:5]
rank = 0
stamp_size = 64
pixel_scale = 0.2
snr_cut_s = float(snr_s)
del_bin = int(del_bin)

ini_path = "%s/work/envs/envs.dat"%my_home
path_items = tool_box.congif(ini_path,['get','get','get','get'], [['selection_bias', "dimmerm_path", '1'],
                                                                  ['selection_bias', "dimmerm_path_result", '1'],
                                                                  ['selection_bias', "dimmerm_path_para", '1'],
                                                                  ['selection_bias', "dimmerm_path_pic", '1']])

total_path = path_items[0]
result_path = path_items[1]
para_path = path_items[2]
pic_path = path_items[3]

shear_input = numpy.load(para_path+"shear.npz")
fg1 = shear_input['arr_0']
fg2 = shear_input['arr_1']

path = result_path + "data/"
for s in range(int(scale)):
    if scale == 1:
        data_cache_path = path + "data_%d.hdf5"%rank
    else:
        data_cache_path = path + 'data_%d_%d.hdf5'%(rank, s)
    f = h5py.File(data_cache_path, 'r')
    if s == 0:
        data = f["/data"].value
    else:
        data = numpy.row_stack((data, f["/data"].value))
    f.close()

fq = Fourier_Quad(stamp_size, 123)


flux = data[:, 7]
detected = flux > 0

# F_Q data
FG1 = data[:, 2]
FG2 = data[:, 3]
FN = data[:, 4]
FU = data[:, 5]
FV = data[:, 6]

prop = lsstetc.ETC(band='r', pixel_scale=pixel_scale, stamp_size=stamp_size, nvisits=180)
noise_sig = prop.sigma_sky

# flux
flux = data[:, 7]/noise_sig
# half_light_flux
hflux = data[:, 8]/noise_sig
# peak
peak = data[:, 9]/noise_sig
# area
area = data[:, 10]
# half_light_area
harea = data[:, 11]
# snr
snr = data[:, 12]
# flux2
flux2 = data[:, 13]
# flux_alt
flux_alt = data[:, 14]

sex_path = total_path + "result/data/sex25_%d_1.5.npz"%rank
sex_data = numpy.load(sex_path)["arr_0"]
cuts_num = 20

mag_auto = sex_data[:, 3]
sex_idx = mag_auto > 0
mag_auto_sort = numpy.sort(-mag_auto[sex_idx])
mag_auto_step = int(len(mag_auto_sort)/cuts_num)
mag_auto_cut = [mag_auto_sort[i*mag_auto_step] for i in range(cuts_num)]

sex_snr = sex_data[:, 0]
sex_snr_sort = numpy.sort(sex_snr[sex_idx])
sex_snr_step = int(len(sex_snr_sort)/cuts_num)
sex_snr_cut = [sex_snr_sort[i*sex_snr_step] for i in range(cuts_num)]

select = {"peak": peak, "flux": flux, "flux2": flux2, "area": area, "snr": snr,
          "hflux": hflux, "harea": harea, "flux_alt": flux_alt, "sex": -mag_auto,"sex_snr": sex_snr}

# detected_label = {"flux": detected, "hflux": detected, "peak": detected, "area": detected, "harea": detected,
#                   "snr": detected, "flux2": detected, "flux_alt": detected, "sex_snr": sex_idx, "sex_area": sex_idx,
#                   "mag_iso": sex_idx, "mag_auto": sex_idx, "mag_petro": sex_idx, "mag_win": sex_idx}
#
res_arr = numpy.zeros((3, 2))
sp = res_arr.shape

ssnr = select[ch]
idxs = ssnr > snr_cut_s

G1 = FG1[idxs]
G2 = FG2[idxs]
N = FN[idxs]
U = FU[idxs]
DE1 = N + U
DE2 = N - U
# weight = ssnr[idxs&idxe]**wei_pow
# if wei_pow == 0:
weight = 1

num = len(G1)

g1_h, g1_h_sig = fq.fmin_g_new(G1, DE1, bin_num=8, ig_num=del_bin)
g2_h, g2_h_sig = fq.fmin_g_new(G2, DE2, bin_num=8, ig_num=del_bin)
print(num, g1_h, g1_h_sig, g2_h, g2_h_sig, fg1[rank], fg2[rank])
numpy.savez(total_path+"%s_PDF.npz"%ch, G2, G2-g2_h*DE2, G2-fg2[rank]*DE2)