import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
from mpi4py import MPI
import numpy
import time


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


dfy_path = "/home/dfy/data/CFHTLens/split-point/"
my_path = "/mnt/perc/hklee/CFHT/catalog/cfht_cata/field_dat/"

pre_fields = tool_box.field_dict("/mnt/perc/hklee/CFHT/catalog/cfht_cata/nname.dat")[1]

my_field = tool_box.allot(pre_fields, cpus)[rank]

nms = ["w1p4p0", "w1p3p0"]

for field_nm in my_field:

    my_cata = numpy.loadtxt(my_path + field_nm + ".dat")

    if field_nm == "w1p4p0":
        dfy_field_nm = "w1p4m0"
    elif field_nm == "w1p3p0":
        dfy_field_nm = "w1p3m0"
    else:
        dfy_field_nm = field_nm
    dfy_cata = numpy.loadtxt(dfy_path + dfy_field_nm + ".dat")

    diff_pos = dfy_cata[:, :2] - my_cata[:, :2]
    max_diff = diff_pos.max()
    print(field_nm, max_diff)
    if max_diff > 0.00001:
        print("Abnormal field %s (%f)"%(field_nm,max_diff))

