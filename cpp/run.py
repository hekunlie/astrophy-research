from subprocess import Popen
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import numpy

num = 50
a = numpy.zeros((num, 1))

sources = ["ptsb"]

for source in sources:
    while True:
        if a.sum() == num:
            break
        for i in range(num):
            if os.path.exists("%s/work/test/job/%s/finish_%d.dat"%(my_home, source, i)):
                a[i, 0] = 1
    print("Run measurement")
    cmd = "mpiexec -n 50 %s/work/cpp/measure_%s_1.5"%(my_home, source)
    run = Popen(cmd, shell=True)
    run.wait()

    cmd = "mpiexec -n 50 %s/work/cpp/measure_%s_2.0"%(my_home, source)
    run = Popen(cmd, shell=True)
    run.wait()

numpy.savetxt("/home/hkli/work/selection_bias/sym_mc_plot/finish.dat",a)
