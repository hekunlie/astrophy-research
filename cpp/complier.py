from subprocess import Popen
from sys import argv

cmd = "mpicxx -I.  %s /home/hkli/work/cpp/FQlib.cpp -lhdf5 -lgsl -lgslcblas -lfftw3l -lfftw3f -lcfitsio  -lm -o %s -std=c++11"\
      %(argv[1],argv[2])
a = Popen(cmd, shell=True)
a.wait()
print("%s >>> %s"%(argv[1], argv[2]))