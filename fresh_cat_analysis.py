from sys import path
path.append('/home/hkli/work/fourier_quad')
import numpy
import os
import tool_box
import time
import shutil
from mpi4py import MPI
import warnings

# to stack the shear catalogs of each exposure into a file

warnings.filterwarnings("error")

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

with open("/home/hkli/work/envs/envs.dat", "r") as f:
    contents = f.readlines()
for path in contents:
    if "cfht_data_path" in path:
        data_path = path.split("=")[1]
    elif "cfht_field_path" in path:
        field_path = path.split("=")[1]

cfht_dict, fields = tool_box.field_dict(data_path + "nname.dat")
field_pool = tool_box.allot(fields, cpus)

for field in field_pool[rank]:
    expos = list(cfht_dict[field].keys())
    f_path = field_path + field + "/"
    for expo in expos:
        i = 0
        for chip in cfht_dict[field][expo]:
            dat_path = data_path + "%s/result/%s_shear.dat"%(field, chip)
            try:
                temp = numpy.loadtxt(dat_path, skiprows=1)
                if i == 0:
                    data = temp
                else:
                    data = numpy.row_stack((data, temp))
                i += 1
            except:
                print("Empty file: %s/%s/%s"%(field,expo,chip))
        expo_path = f_path + "%s.npz"%expo
        numpy.savez(expo_path, data)


