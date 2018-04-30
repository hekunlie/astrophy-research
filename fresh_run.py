from subprocess import Popen
import time
import os
import shutil


# specify the number of cores used
ncpus = 64

# specify it to where the para.inc is stored
para = "/home/hkli/work/fresh/src_new/para.inc"

# specify it to where the data are stored (don't forget the last "/")
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/w1234/original/"
source_list = data_path + "source.list"

# run the three processes
# one can specify the PROCESS_stage like [1,3]...
PROCESS_stage = [1, 2, 3]

# read the para.inc
f = open(para, "r")
contents = f.readlines()
f.close()

dirs = os.listdir(data_path)
initial_target = ["/astrometry", "/stamps", "/result"]

t = [0, 0, 0]

for i in PROCESS_stage:
    t1 = time.time()

    # initialize the directories, "../astrometry" or " ../result" or "../stamps" (depending on the PROCESS_stage)
    # it will take a while...
    for dir in dirs:
        if "w" in dir:
            initial_path = data_path + dir + initial_target[i-1]
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
