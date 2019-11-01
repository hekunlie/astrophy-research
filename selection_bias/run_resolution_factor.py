from subprocess import Popen
import time
from sys import argv


sources = ["galsim_bright", "galsim_dimmer", "pts_bright", "pts_dimmer"]
filters = ["sex2_1.5", "sex2_2", "sex2_4", "sex4_1.5", "sex4_2", "sex4_4"]

for ss in sources:
    for ff in filters:

        print("----------------------------------------------------------")
        t1 = time.time()
        cmd = "mpirun -n 25 ./resolution %s %s 1000"%(ss,ff)
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print(cmd + " %.2f sec"%(t2-t1))
        print("----------------------------------------------------------\n")