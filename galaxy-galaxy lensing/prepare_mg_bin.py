import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
import h5py
from mpi4py import MPI
from sys import argv
import matplotlib.pyplot as plt
import numpy
import time



areas = [1,3,4]


envs_path = "%s/work/envs/envs.dat"%my_home

gets_item = [["cfht", "cfht_path_catalog", "0"], ["gg_lensing", "ggl_path_data", "0"]]
path_items = tool_box.config(envs_path, ["get", "get"], gets_item)
cata_path, data_path = path_items

h5f_path = data_path + "cata_result_ext_grid.hdf5"

h5f = h5py.File(h5f_path)

for area_id in areas:
    # the G*\Simg_crit
    h5f_mg = h5py.File(data_path + "w_%d_/10.hdf5","r")
    mg_crit = h5f_mg[""]





