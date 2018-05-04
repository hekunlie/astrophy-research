from subprocess import Popen
import time
import os
import shutil
from sys import argv


# specify the number of cores used
ncpus = 64

# specify it to where the para.inc is stored
para = "/home/hkli/work/fresh/src_new/para.inc"

# specify it to where the data are stored (don't forget the last "/")
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/w1234/original/"
source_list = data_path + "source.list"

# run the three processes
# one can specify the PROCESS_stage like [1,3], which comes from the commands from command inputs
PROCESS_stage = []
for i in range(len(argv)):
    if 0 < i < 4:
        c = int(float(argv[i]))
        if c in [1, 2, 3]:
            PROCESS_stage.append(c)
        else:
            print("Wrong PROCESS_stage input: %s which should be one of [1, 2, 3]"%argv[i])
            exit()

# read the para.inc
f = open(para, "r")
contents = f.readlines()
f.close()

files = os.listdir(data_path)
initial_target = ["/astrometry", "/stamps", "/result"]

t = [0, 0, 0]

for i in PROCESS_stage:
    t1 = time.time()

    # initialize the directories, "../astrometry" or " ../result" or "../stamps" (depending on the PROCESS_stage)
    # it will take a while...
    for file in files:
        if "w" in dir:
            initial_path = data_path + file + initial_target[i-1]
            print("Deleting the %s"%initial_path)
            shutil.rmtree(initial_path, ignore_errors=True)
            os.mkdir(initial_path)

    # change the parameter "PROCESS_stage" in the "para.inc"
    # the blank space before "parameter" is required by fortran77
    contents[10] = "	parameter (PROCESS_stage=%d)\n"%i
    f = open(para, "w")
    f.writelines(contents)
    f.close()

    # compiling
    cmd = "mpif77 *.f -o main -mcmodel=medium -lcfitsio"
    a = Popen(cmd, shell=True)
    a.wait()
    # run code
    cmd = "mpiexec -n %d ./main "%ncpus + source_list
    a = Popen(cmd, shell=True)
    a.wait()

    t2 = time.time()
    t[i-1] = (t2-t1)/3600

print("PROCESS 1: %.2f H, PROCESS 2: %.2f H, PROCESS 3: %.2f H"%(t[0], t[1], t[2]))
