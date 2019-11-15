from subprocess import Popen
from sys import argv
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]

compiler = "mpicxx"
cpp_std = "-std=c++11"

lib_num = len(argv)
# code name
code_name = argv[1]
# program name after compiling
final_program_name = argv[2]

# the libs will be included
# "FQlib", "hk_iolib" are available now
link_libs = [argv[i] for i in range(2, lib_num)]

# the external libs for "hk_iolib"
link_hk_iolib = "-lhdf5 -lcfitsio "
# the external libs for "FQlib"
link_FQlib = "-lgsl -lgslcblas -lfftw3l -lfftw3f "

links = ""
cmd = "%s -I.  %s "%(compiler,code_name)
for i in range(2,lib_num):
      if "hk_iolib" in argv[i]:
            links += link_hk_iolib
            cmd += "%s/work/mylib/%s.cpp "%(my_home, argv[i])
      if "FQlib" in argv[i]:
            links += link_FQlib
            cmd += "%s/work/mylib/%s.cpp " %(my_home, argv[i])

cmd += " %s -lm -o %s %s "%(links,final_program_name,cpp_std)
print(cmd)
a = Popen(cmd, shell=True)
a.wait()
print("%s >>> %s"%(argv[1], argv[2]))