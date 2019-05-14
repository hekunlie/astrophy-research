from subprocess import Popen
import time
from sys import argv
import os

foreground = argv[2]
area = [1, 2, 3, 4]
radius_num = 13
if argv[1] == "search":
    for area_id in area:
        for radius_id in range(radius_num):
            cmd = "mpirun -np 1 ./ggl_test %d %d %s"%(area_id, radius_id, foreground)
            t1 = time.time()
            a = Popen(cmd, shell=True)
            a.wait()
            t2 = time.time()
            print("Finish %.2f sec"%(t2-t1))

    print("\nMeasure\n")
    for area_id in area:
        cmd = "mpirun -np %d python ggl_sym.py %d %s"%(radius_num,area_id, foreground)
        t1 = time.time()
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print("Finish %.2f sec"%(t2-t1))

elif argv[1] == "measure":
    print("\nMeasure\n")
    for area_id in area:
        cmd = "mpirun -np %d python ggl_sym.py %d %s"%(radius_num,area_id, foreground)
        t1 = time.time()
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print("Finish %.2f sec"%(t2-t1))