from subprocess import Popen
import time
from sys import argv


foreground = argv[2]
area = [1, 2, 3, 4]
radius_num = 13
cpus = 40
if argv[1] == "search":
    for area_id in area:
        for radius_id in range(radius_num):
            cmd = "mpirun -np 1 ./ggl_test %d %d %s"%(area_id, radius_id, foreground)
            t1 = time.time()
            a = Popen(cmd, shell=True)
            a.wait()
            t2 = time.time()
            print("Finish %.2f sec"%(t2-t1))

elif argv[1] == "measure_all":
    print("\nMeasure\n")
    for area_id in area:
        cmd = "mpirun -np %d python ggl_sym.py %d %s"%(radius_num,area_id, foreground)
        t1 = time.time()
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print("Finish %.2f sec"%(t2-t1))

elif argv[1] == "gamma":
    print("\nMeasure\n")
    for radius_id in range(radius_num):
        cmd = "mpirun -n %d  ./ggl_sym 1 %s %s 1 2 3 4"%(cpus, foreground, radius_id)
        t1 = time.time()
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print("Finish %.2f sec"%(t2-t1))

elif argv[1] == "crit":
    print("\nMeasure\n")
    for radius_id in range(radius_num):
        cmd = "mpirun -n %d ./ggl_sym 2 %s %s 1 2 3 4"%(cpus, foreground, radius_id)
        t1 = time.time()
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print("Finish %.2f sec"%(t2-t1))