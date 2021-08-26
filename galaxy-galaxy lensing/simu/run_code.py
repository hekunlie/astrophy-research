from subprocess import Popen
import os
import shutil
import time


data_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_1"

cmd = "mpirun -np 20 python %s/NFW_shear_map_continue_source_plane.py 0 3000000 0"%data_path
a = Popen(cmd, shell=True)
a.wait()

cmd = "mpirun -np 50 python %s/NFW_shear_map_continue_source_plane.py 0 3000000 1"%data_path
a = Popen(cmd, shell=True)
a.wait()

for i in range(50):
    cmd = "mpirun -np 50 ./simu /home/hklee/work/Galaxy_Galaxy_lensing_test %d 300 sheared" %i
    a = Popen(cmd, shell=True)
    a.wait()


for i in range(20):
    cmd = "mpirun -np 50 ./simu /home/hklee/work/Galaxy_Galaxy_lensing_test %d 300 non_sheared" %i
    a = Popen(cmd, shell=True)
    a.wait()
