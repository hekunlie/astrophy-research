from subprocess import Popen
import time
from sys import argv


fore_source_nm = argv[1]
area = [int(argv[i]) for i in range(2, len(argv))]
radius_num = 13
print(area)
time.sleep(2)

cpus = [2 + i*3 for i in range(radius_num)]

for area_id in area:
    for radius_id in range(radius_num):
        t1 = time.time()
        cmd = "mpirun -n %d ./ggl_cal %d %d %s"%(cpus[radius_id], area_id, radius_id, fore_source_nm)
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print(cmd + " %.2f sec"%(t2-t1))
        print("\n\n\n")