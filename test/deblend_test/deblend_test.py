import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
import c4py
import time

expo_name = argv[1]
expo_path = "/mnt/perc/hklee/CFHT/CFHT_cat_4_20_2021/cat_hdf5/%s"%expo_name

sep_pix_thresh = 12
sep_z_thresh = 0.2


h5f = h5py.File(expo_path,"r")
expo_data = h5f["/data"][()]
h5f.close()


expo_ichip_all = expo_data[:,16].astype(dtype=numpy.intc)

expo_xc_all = expo_data[:, 18]
expo_yc_all = expo_data[:, 19]
expo_z_all = expo_data[:, 10]
expo_radius_all = numpy.sqrt(expo_data[:,25]/numpy.pi)

print(expo_ichip_all.min(), expo_ichip_all.max())

nums_all = numpy.zeros((expo_ichip_all.max()+1, ), dtype=numpy.intc)
nums = numpy.zeros((expo_ichip_all.max()+1, ), dtype=numpy.intc)
nums_i = numpy.zeros((expo_ichip_all.max()+1, ), dtype=numpy.intc)
nums_i_all = numpy.zeros((expo_ichip_all.max()+1, ), dtype=numpy.intc)


for i in range(expo_ichip_all.max()+1):
    if i == 0:
        idx_ic = expo_ichip_all > -1
    else:
        idx_ic = expo_ichip_all == i

    expo_ichip = expo_ichip_all[idx_ic]
    expo_xc = expo_xc_all[idx_ic]
    expo_yc = expo_yc_all[idx_ic]
    expo_z = expo_z_all[idx_ic]
    expo_radius = expo_radius_all[idx_ic]

    cata_len_1 = expo_xc.shape[0]

    gal_label = numpy.arange(0, cata_len_1)


    t1 = time.time()
    labels = c4py.deblend(expo_xc, expo_yc, expo_z, expo_radius, expo_ichip, 4, sep_z_thresh)
    t2 = time.time()
    idx1 = labels < 0
    print(idx1.sum(),cata_len_1)
    # print(gal_label[idx1])
    print("%.2f sec"%(t2-t1))

    nums[i] = idx1.sum()
    nums_all[i] = cata_len_1


    labels_i = c4py.deblend_i(expo_xc, expo_yc, expo_z, expo_ichip, sep_pix_thresh, sep_z_thresh)

    t3 = time.time()
    idx1_i = labels < 0
    print(idx1_i.sum(),cata_len_1)
    # print(gal_label[idx1_i])
    print("%.2f sec"%(t3-t2))
    nums_i[i] = idx1_i.sum()
    nums_i_all[i] = cata_len_1

print(nums)
print(nums_all)
print(nums[1:].sum(),nums[1:].sum()/nums_all[0])

print(nums_i)
print(nums_i_all)
print(nums_i[1:].sum(),nums_i[1:].sum()/nums_all[0])
# labels_test = numpy.ones((cata_len_1,), dtype=numpy.intc)
# for i in range(cata_len_1):
#
#     radius_i = expo_radius[i]
#     xc_i, yc_i = expo_xc[i], expo_yc[i]
#     z_i = expo_z[i]
#
#     for j in range(i, cata_len_1):
#
#         r2 = radius_i + expo_radius[j]
#         dx = numpy.abs(xc_i - expo_xc[j]) - r2
#
#         if dx <= 0:
#             dy = numpy.abs(yc_i - expo_yc[j]) - r2
#             dz = numpy.abs(z_i - expo_z[j])
#             if dy <= 0 and dz >= sep_z_thresh:
#                 labels_test[i] = -1
#                 labels_test[j] = -1
#                 # print(i, j)
#                 # print(expo_xc[i], expo_xc[j], numpy.abs(xc_i - expo_xc[j]))
#                 # print(expo_yc[i], expo_yc[j], numpy.abs(yc_i - expo_yc[j]))
#                 # print(r2)
#                 # print(expo_z[i], expo_z[j], dz)
#                 # print("\n")
#
# t4 = time.time()
# idx2 = labels_test < 0
# print(idx2.sum(), idx2.shape[0])
# print("%.2f sec" % (t4 - t3))
#
#
# idx = idx1 & idx2
# print(idx.sum())
# idx = idx1 | idx2
# print(idx.sum())
# print(gal_label[idx1])
# print(gal_label[idx2])
# print(gal_label[idx1]-gal_label[idx2])



