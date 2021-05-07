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
pz_path = "/mnt/perc/hklee/CFHT/catalog/CFHTLens_Pz.hdf5"


h5f = h5py.File(expo_path,"r")
expo_data = h5f["/data"][()]
h5f.close()
expo_ra = expo_data[:,0]
expo_dec = expo_data[:,1]
expo_z = expo_data[:,10]
cata_len_1 = expo_ra.shape[0]

expo_ra_min, expo_ra_max = expo_ra.min(), expo_ra.max()
dra = expo_ra_max - expo_ra_min

expo_dec_min, expo_dec_max = expo_dec.min(), expo_dec.max()
ddec = expo_dec_max - expo_dec_min


h5f = h5py.File(pz_path,"r")
pz_data = h5f["/data"][()]
h5f.close()

idx1 = pz_data[:,0] >= expo_ra_min - dra*0.05
idx2 = pz_data[:,0] <= expo_ra_max + dra*0.05
idx3 = pz_data[:,1] >= expo_dec_min - ddec*0.05
idx4 = pz_data[:,1] <= expo_dec_max + ddec*0.05
idx = idx1 & idx2 & idx3 & idx4

pz_ra = pz_data[:,0][idx]
pz_dec = pz_data[:,1][idx]
pz_z = pz_data[:,2][idx]
cata_len_2 = pz_ra.shape[0]

t1 = time.time()

seq_label, check_label = c4py.cata_match(expo_ra, expo_dec, pz_ra, pz_dec, 0.0003)

t2 = time.time()

print("%.2f sec"%(t2-t1))

pair_test = numpy.zeros((cata_len_1, 3))

for i in range(cata_len_1):
    tag = seq_label[i]
    pair_test[i] = pz_ra[tag], pz_dec[tag], pz_z[tag]
    if i < 10:
        print(expo_ra[i], pz_ra[tag])
diff_ra = expo_ra - pair_test[:,0]
diff_dec = expo_dec - pair_test[:,1]


print(diff_ra.min(), diff_ra.max())
print(diff_dec.min(), diff_dec.max())


pair_test_n = numpy.zeros((cata_len_1, 3))
t3 = time.time()
for i in range(cata_len_1):
    diff_rad = numpy.abs(pz_ra - expo_ra[i]) + numpy.abs(pz_dec - expo_dec[i])

    idx = diff_rad == diff_rad.min()

    pair_test_n[i] = pz_ra[idx],pz_dec[idx],pz_z[idx]

t4 = time.time()
print("%.2f sec"%(t4-t3))

diff_ra = expo_ra - pair_test_n[:,0]
diff_dec = expo_dec - pair_test_n[:,1]


print(diff_ra.min(), diff_ra.max())
print(diff_dec.min(), diff_dec.max())


