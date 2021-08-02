from subprocess import Popen
import os
import shutil
import time

for i in range(50):
    cmd = "mpirun -np 10 ./simu /home/hklee/work/Galaxy_Galaxy_lensing_test %d 10 sheared" %i
    a = Popen(cmd, shell=True)
    a.wait()
