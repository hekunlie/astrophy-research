import matplotlib
matplotlib.use('Agg')
from sys import path, argv
path.append('/home/hkli/work/fourier_quad')
from Fourier_Quad import Fourier_Quad
import numpy
import time
import os
import matplotlib.pyplot as plt
import tool_box
import copy


data_name, g1num, g2num, bin_num, thresh = argv[1],argv[2],argv[3],argv[4], float(argv[5])
g1num, g2num, bin_num = int(g1num), int(g2num), int(bin_num)

with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_res" in path:
        result_path = path.split("=")[1]
field_path = result_path + "field/"

filter_path = result_path + "field/filtered.dat"
filter_exist = os.path.exists(filter_path)

fq = Fourier_Quad(48,123)
cache_path = result_path + data_name
print(cache_path)
arr = numpy.load(cache_path)["arr_1"]
est_g1 = arr[:g1num, 0]
fd_g1 = arr[:g1num, 6]
est_g2 = arr[:g2num, 3]
fd_g2 = arr[:g2num, 7]

dg1 = fd_g1[1] - fd_g1[0]
dg2 = fd_g2[1] - fd_g2[0]

files = []
if filter_exist:
    with open(filter_path, "r") as f:
        contents = f.readlines()

    for c in contents:
        name = c.split("\n")[0] + "_shear.dat.npz"
        files.append(name)
    print("Filtered data: %d" % len(files))
else:
    file_list = os.listdir(field_path)
    count = 0
    for i in range(len(file_list)):
        if "shear.dat.npz" in file_list[i]:
            files.append(file_list[i])
            count += 1
    print("Unfiltered data: %d"%count)

total_chi_1 = []
total_chi_2 = []
field_name = []

for name in files:
    if "w" in name:
        chi_1s = []
        chi_2s = []
        new_path = field_path + name
        field_data = numpy.load(new_path)["arr_0"]
        n_star = field_data[:, 3]
        idx = n_star >= 16
        fsnr = field_data[:, 10][idx]
        fg1 = field_data[:, 14][idx]
        fg2 = field_data[:, 15][idx]
        FG1 = field_data[:, 16][idx]
        FG2 = field_data[:, 17][idx]
        FN = field_data[:, 18][idx]
        FU = field_data[:, 19][idx]

        idx_thres = fsnr >= thresh

        for i in range(g1num):
            idx11 = fg1 >= fd_g1[i] - dg1 / 2
            idx12 = fg1 <= fd_g1[i] + dg1 / 2
            mg1 = FG1[idx11&idx12&idx_thres]
            mn = FN[idx11&idx12&idx_thres]
            mu = FU[idx11&idx12&idx_thres]
            mg1_bins = tool_box.set_bin(mg1, bin_num)
            chi_1 = fq.G_bin(mg1, mn, mu, est_g1[i], 1, mg1_bins, 0)
            chi_1s.append(chi_1)

        for j in range(g2num):
            idx21 = fg2 >= fd_g2[j] - dg2 / 2
            idx22 = fg2 <= fd_g2[j] + dg2 / 2
            mg2 = FG2[idx21&idx22&idx_thres]
            mn = FN[idx21&idx22&idx_thres]
            mu = FU[idx21&idx22&idx_thres]
            mg2_bins = tool_box.set_bin(mg2, bin_num)
            chi_2 = fq.G_bin(mg2, mn, mu, est_g2[j], 2, mg2_bins, 0)
            chi_2s.append(chi_2)
        total_chi_1.append(chi_1s)
        total_chi_2.append(chi_2s)
        field_name.append(name)


thresh1 = 3
thresh2 = 2.5
over_sig_per = [0.3, 0.3]
fg_bin_num = [g1num, g2num]
total_chis = [numpy.array(total_chi_1), numpy.array(total_chi_2)]
field_tag = numpy.ones((len(field_name), 1))
masks = copy.deepcopy(total_chis)
numpy.savez("/mnt/ddnfs/data_users/hkli/CFHT/result/field/chi.npz",total_chis[0],total_chis[1])


for i in range(2):
    lb = numpy.arange(0, len(field_name))
    # to exclude the outlier of the total chi square
    total_chi_c = total_chis[i].copy()
    mask = numpy.ones_like(total_chi_c)
    sig = numpy.std(total_chi_c)
    idx = total_chi_c < thresh1*sig
    mask[idx] = 0
    overflow = numpy.sum(mask, axis=1)
    idx = overflow > over_sig_per[0]*fg_bin_num[i]
    # label the excluded fields
    field_tag[lb[idx]] = 0
    masks[i][lb[idx]] = -5

    # to exclude the outliers bin-wise
    idx = overflow <= over_sig_per[0]*fg_bin_num[i]
    field_chi = numpy.sum(total_chi_c[lb[idx]], axis=1)[:,numpy.newaxis]
    field_sig = numpy.std(field_chi)
    scale = (field_chi - numpy.mean(field_chi))/field_sig
    idxs = scale > thresh2
    lb_c = lb[idx][:,numpy.newaxis]
    field_tag[lb_c[idxs]] = 0
    masks[i][lb_c[idxs]] = -15

filtered = []
excluded = []
for i in range(len(field_name)):
    name = field_name[i].split("_")[0] + "\n"
    if field_tag[i, 0] != 0:
        filtered.append(name)
    else:
        excluded.append(name)
        print(i, field_name[i])

filtered_path = field_path + "filtered.dat"
with open(filtered_path, "w") as f:
    f.writelines(filtered)

excluded_path = field_path + "excluded.dat"
with open(excluded_path, "w") as f:
    f.writelines(excluded)

plt.figure(figsize=(6,10))
plt.subplot(121)
plt.imshow(total_chis[0])
plt.colorbar()
plt.subplot(122)
plt.imshow(total_chis[1])
plt.colorbar()
png = field_path + "chi.png"
plt.savefig(png)
plt.close()

plt.figure(figsize=(6,10))
plt.subplot(121)
plt.imshow(masks[0])
plt.colorbar()
plt.subplot(122)
plt.imshow(masks[1])
plt.colorbar()
png = field_path + "chi_mask.png"
plt.savefig(png)
plt.close()



