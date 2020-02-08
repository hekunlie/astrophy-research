from subprocess import Popen
from sys import argv

data_path = argv[1]
data_type = argv[2]
data_num = int(argv[3])

data_row = [4,5,6]

# non-weighted
cmd = "mpirun -np 15 ./calculate %s %s 15 %d 0 0" % (data_path, data_type, data_num)
a = Popen(cmd, shell=True)
a.wait()

# weighted
for i in data_row:
    for j in range(1,5):
        cmd = "mpirun -np 15 ./calculate %s %s 15 %d %d %d"%(data_path, data_type,data_num, i, j)
        a = Popen(cmd, shell=True)
        a.wait()


