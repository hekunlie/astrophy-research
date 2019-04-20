from subprocess import Popen
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import numpy

num = 50

sources = ["dimmerm"]

jobs = numpy.zeros((num, 1))
for source in sources:

    while True:
        for i in range(50):
            if os.path.exists("/home/hkli/work/test/job/%s/finish_%d.dat"%(source,i)):
                jobs[i] = 1
        if jobs.sum() == num:
            break

    print("Run measurement")
    cmd = "mpiexec -n 40 %s/work/cpp/measure_%s_1.5"%(my_home, source)
    run = Popen(cmd, shell=True)
    run.wait()

    cmd = "mpiexec -n 40 %s/work/cpp/measure_%s_2.0"%(my_home, source)
    run = Popen(cmd, shell=True)
    run.wait()

    cmd = "mpiexec -n 40 %s/work/cpp/measure_%s_4.0"%(my_home, source)
    run = Popen(cmd, shell=True)
    run.wait()

    with open("/home/hkli/work/test/job/%s/finish.dat"%source) as f:
        f.write("f")
