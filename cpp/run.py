from subprocess import Popen
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import numpy

num = 50

sources = ["pts"]

for source in sources:
    while True:
        if os.path.exists("/home/hkli/work/cpp/finish.dat"):
            break
    print("Run measurement")
    cmd = "mpiexec -n 50 %s/work/cpp/measure_%s_1.5"%(my_home, source)
    run = Popen(cmd, shell=True)
    run.wait()

    cmd = "mpiexec -n 50 %s/work/cpp/measure_%s_2.0"%(my_home, source)
    run = Popen(cmd, shell=True)
    run.wait()

with open("/home/hkli/work/selection_bias/sym_mc_plot/finish.dat") as f:
    f.write("f")
