from subprocess import Popen
from sys import argv
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]

cmd = "mpicxx -I.  %s %s/work/cpp/mylib/FQlib.cpp -lhdf5 -lgsl -lgslcblas -lfftw3l -lfftw3f -lcfitsio  -lm -o %s -std=c++11"\
      %(argv[1],my_home, argv[2])
a = Popen(cmd, shell=True)
a.wait()
print("%s >>> %s"%(argv[1], argv[2]))