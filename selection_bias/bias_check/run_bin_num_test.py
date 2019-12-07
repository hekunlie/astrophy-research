from subprocess import Popen
from sys import argv

bin_num = [2,4,8,12,16,20,40,60,80,120,160,200]

src_path = argv[1]
numprocs = 17
shear_num = 17
data_type = "noisy"
src_num = 10000

for bn in bin_num:
    for j in range(9):
        cmd = "mpirun -np %d ./calculate %s %d %s %d %d %d"%(numprocs, src_path, shear_num, data_type, src_num, bn, j)
        print(cmd)
        a = Popen(cmd, shell=True)
        a.wait()
