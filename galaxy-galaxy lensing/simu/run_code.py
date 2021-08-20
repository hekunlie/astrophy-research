from subprocess import Popen
import os
import shutil
import time



cmd = "mpirun -np 50 python /home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_1/NFW_shear_map_continue_source_plane.py 5000000 0"
a = Popen(cmd, shell=True)
a.wait()
cmd = "mpirun -np 50 python /home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_1/NFW_shear_map_continue_source_plane.py 5000000 1"
a = Popen(cmd, shell=True)
a.wait()

for i in range(50):
    cmd = "mpirun -np 50 ./simu /home/hklee/work/Galaxy_Galaxy_lensing_test %d 500 sheared" %i
    a = Popen(cmd, shell=True)
    a.wait()
for i in range(50):
    cmd = "mpirun -np 50 ./simu /home/hklee/work/Galaxy_Galaxy_lensing_test %d 500 non_sheared" %i
    a = Popen(cmd, shell=True)
    a.wait()
