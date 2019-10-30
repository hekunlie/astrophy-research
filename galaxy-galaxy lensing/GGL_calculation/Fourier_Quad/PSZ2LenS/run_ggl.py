from subprocess import Popen
import time
from sys import argv


fore_source_nm = argv[1]
radius_num = int(argv[2])
area = [int(argv[i]) for i in range(3, len(argv))]
print(area)
time.sleep(2)

cpus = [1 for i in range(radius_num)]

for area_id in area:
    for radius_id in range(radius_num):

        print("----------------------------------------------------------")
        t1 = time.time()
        cmd = "mpirun -n %d ./ggl_cal %d %d %s"%(cpus[radius_id], area_id, radius_id, fore_source_nm)
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print(cmd + " %.2f sec"%(t2-t1))
        print("----------------------------------------------------------\n")
