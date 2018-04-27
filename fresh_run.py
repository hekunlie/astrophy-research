from subprocess import Popen
import time
import os
import shutil


code_path = "/home/hkli/work/fresh/src_new/"
para = code_path + "para.inc"

data_path = "/mnt/ddnfs/data_users/hkli/CFHT/w1234/original/"
source_list = data_path + "source.list"

f = open(para, "r")
contents = f.readlines()
f.close()

# initializing the directories "../astrometry, ../result, ../stamps"
# it will take a while...
dirs = os.listdir(data_path)
for dir in dirs:
    if "w" in dir:
        astro = data_path + dir + "/astrometry"
        shutil.rmtree(astro, ignore_errors=True)
        os.mkdir(astro)

        result = data_path + dir + "/result"
        shutil.rmtree(result, ignore_errors=True)
        os.mkdir(result)

        stamps = data_path + dir + "/stamps"
        shutil.rmtree(stamps, ignore_errors=True)
        os.mkdir(stamps)

# run the three processes
PROCESS_stage = [1, 2, 3]
t = []
for i in PROCESS_stage:
    t1 = time.clock()

    # change the PROCESS_stage
    contents[10] = "	parameter (PROCESS_stage=%d)\n"%i
    f = open(para, "w")
    f.writelines(contents)
    f.close()

    # compiling
    cmd = "mpif77 *.f -o main -mcmodel=medium -lcfitsio"
    a = Popen(cmd, shell=True)
    a.wait()
    # run code with 64 threads
    cmd = "mpiexec -n 64 ./main " + source_list
    a = Popen(cmd, shell=True)
    a.wait()

    t2 = time.clock()
    t.append(t2-t1)

print(t)
