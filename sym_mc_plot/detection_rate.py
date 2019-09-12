import numpy
import h5py
import os

fs = os.listdir("../paper_data/")

files = [i for i in fs if "bright" in i or "dimmer" in i]

for nm in files:
    rate_det = []
    for i in range(5):
        mask_path = "../paper_data/%s/result/data/sex2_1.5/mask_%d.hdf5"%(nm,i)
        h5f = h5py.File(mask_path, "r")
        mask = h5f["/data"].value
        idx = mask == 1
        rate_det.append(idx.sum()/mask.shape[0])
    print(nm, rate_det)