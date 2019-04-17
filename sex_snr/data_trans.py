import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append("%s/work/mylib/"%my_home)
import tool_box
from subprocess import Popen
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

ts = time.time()

data_path = "/mw/w1234/original/"
nm_path = data_path + "nname.dat"
fields = tool_box.field_dict(nm_path)[1]

missions = tool_box.allot(fields, cpus)[rank]

store_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/"

for dirs in ["result", "result_int", "result_ext"]:
    for field in missions:
        try:
            cmd = "scp -r /mw/w1234/original/%s/%s/ " \
                  "hkli@202.120.32.231:/mnt/ddnfs/data_users/hkli/CFHT/catalog/%s/"\
                  %(field, dirs, field)
            # print(cmd)
            a = Popen(cmd, shell=True)
            a.wait()
        except:
            print("FAILED %s"%field)

te = time.time()
print(te-ts)
