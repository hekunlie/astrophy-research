from subprocess import Popen
from sys import argv

cpus = int(argv[1])
rotation_lb = [int(argv[i]) for i in range(2,len(argv))]

for i in rotation_lb:
    cmd = "mpirun -np %d ./noise_test %d"%(cpus,i)
    a = Popen(cmd, shell=True)
    a.wait()
