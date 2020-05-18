from subprocess import Popen
from sys import argv

data_type = argv[1]
data_row = [0,4,5,6]
for i in data_row:
    cmd = "mpirun -np 15 ./calculate /mnt/perc/hklee/bias_check/weight_test/result %s 15 %d"%(data_type, i)
    a = Popen(cmd, shell=True)
    a.wait()
