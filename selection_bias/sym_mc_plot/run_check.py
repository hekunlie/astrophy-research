from subprocess import Popen
from sys import argv
st, ed = int(argv[1]), int(argv[2])

for i in range(st,ed):
    cmd = "mpirun -np 1 python prepare_data.py /coma/hklee %d" %i
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "mpirun -np 15 python cutoff_cfht.py /coma/hklee/selection_bias %d" %i
    a = Popen(cmd, shell=True)
    a.wait()
