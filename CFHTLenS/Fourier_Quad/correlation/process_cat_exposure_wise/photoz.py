import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
import h5py
from tqdm import tqdm

# prepare the photoz cat from the full Pz
# ra, dec, Z_B, Z_B_MIN, Z_B_MAX, ODDS, P(Z)(70 columns)
cata_path = "/mnt/perc/hklee/CFHT/catalog"

with open(cata_path + "/CFHTLens_Pz.tsv", "r") as f:
    cc = f.readlines()
row = len(cc) - 1
col = 76
print(row, " obj")

zbin = numpy.array([0.025 + i * 0.05 for i in range(70)])

if argv[1] == "prepare":

    data_1 = numpy.zeros((row, 8), dtype=numpy.float32)
    data_2 = numpy.zeros((row, 70), dtype=numpy.float32)

    pz_ab = []
    data_ab = []
    for i in tqdm(range(row)):

        cc_s = cc[i+1].replace("\"","").replace(" ","").split(",")

        if len(cc_s) != col:
            print("Line %d. Column does not match!!!"%i)
            print(len(cc_s))
            print(cc[i+1])
            print(cc_s)
            exit()

        for tag, element in enumerate(cc_s[:6]):
            data_1[i, tag] = float(element)
        for tag, element in enumerate(cc_s[6:]):
            data_2[i, tag] = float(element)

        # the expectation of Z
        # some of the P(z) does not be normalized to 1
        pz = data_2[i]
        # some of the Sum may equal zero!
        pz_sum = pz.sum()
        # if pz_sum < 0.9:
        #     print("%d P(z) wrong!!!"%i)
        #     print(pz,pz_sum, data_1[i,3])
        #     pz_ab.append(pz)
        #     data_ab.append(data_1[i])

        pz_norm = pz / pz_sum
        data_1[i, 6] = numpy.sum(pz_norm*zbin)
        data_1[i, 7] = pz_sum

    # pz_ab = numpy.array(pz_ab)
    # data_ab = numpy.array(data_ab)
    # numpy.savez("cache.npz", pz_ab, data_ab)


    h5f = h5py.File(cata_path + "/CFHT_pz.hdf5", "w")
    h5f["/data"] = data_1
    h5f["/data"].attrs["/col_name"] = ["ra, dec, Z_B, Z_B_MIN, Z_B_MAX, ODDS, Z_expect, Pz_sum"]
    h5f["/pz"] = data_2
    h5f.close()

else:
    # test
    h5f = h5py.File(cata_path + "/CFHT_pz.hdf5", "r")
    data_1 = h5f["/data"][()]
    data_2 = h5f["/pz"][()]
    h5f.close()

    # row_label = numpy.random.randint(0, row, 20000)
    # temp = numpy.zeros((col,), dtype=numpy.float32)
    # for i in row_label:
    #     cc_s = cc[i + 1].replace(",", "\t").split("\t")
    #     if len(cc_s) != col:
    #         print("Line %d. Column does not match!!!" % i)
    #         print(len(cc_s))
    #         print(cc[i + 1])
    #         print(cc_s)
    #         exit()
    #     for tag, element in enumerate(cc_s):
    #         temp[tag] = float(element)
    #
    #     diff = temp - data_1[i]
    #     if diff.max() >= 0.0001 or diff.min() <= -0.0001:
    #         print("%d line wrong"%i)
    #         print(data[i],"\n")
    #         print(cc_s,"\n")
    #         print(temp,"\n")
    #
    # # check the last line
    # cc_s = cc[-1].replace(",", "\t").split("\t")
    # for tag, element in enumerate(cc_s):
    #     temp[tag] = float(element)
    #
    # diff = temp - data[-1]
    # if diff.max() >= 0.0001 or diff.min() <= -0.0001:
    #     print("Last line wrong")
    #     print(data[-1], "\n")
    #     print(cc_s, "\n")
    #     print(temp, "\n")




