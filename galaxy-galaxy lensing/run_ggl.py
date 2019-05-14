from subprocess import Popen
import time
from sys import argv
import numpy

cpus = int(argv[1])
area = [1,3,4]
radius_label = [i for i in range(13)]
radius_start = numpy.linspace(-110, -15, 13)
radius_end = -radius_start
point_num = [int(3*i) for i in radius_end]
print(radius_start)
print(radius_end)
print(point_num)
for area_id in area:
    for tag, radius_id in enumerate(radius_label):
        t1 = time.time()
        cmd = "mpirun -n %d ./ggl_cal %d %d %f %f %d"%(cpus, area_id, radius_id,
                                                       radius_start[tag], radius_end[tag], point_num[tag])
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print(cmd + " %.2f sec"%(t2-t1))
        print("\n\n\n")