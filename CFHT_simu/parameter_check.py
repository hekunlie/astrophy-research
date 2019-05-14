import platform
if platform.system() == 'Linux':
    import matplotlib
    matplotlib.use('Agg')
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/'%my_home)
import numpy
import matplotlib.pyplot as plt
import tool_box
from mpi4py import MPI
import h5py
import time

sect, source,loops = argv[1], argv[2], int(argv[3])

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

envs_path = "%s/work/envs/envs.dat"%my_home
para_path = tool_box.config(envs_path,["get"],[["%s"%sect,"%s_path_para"%source,"0"]])[0]

para_ini_path = para_path+"para.ini"
paras = tool_box.config(para_ini_path,["get",'get',"get",'get',"get",'get'],
                        [["para","total_num","0"],["para","stamp_size","0"],
                        ["para", "mag_s", "0"],["para","mag_e","0"],
                        ["para", "radius_s", "0"],["para","radius_e","0"]])


pic = para_path + "/pic/ellip_%d.png"%rank
plt.figure(figsize=(16,16))

h5_path = para_path+'para_%d.hdf5'%rank
f = h5py.File(h5_path,"r")


