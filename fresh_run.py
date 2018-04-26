from subprocess import Popen
import time


code_path = "/home/hkli/work/fresh/src_new/"
para_mode = code_path + "para_mode.inc"
para_process = code_path + "para.inc"
source_list = "/mnt/ddnfs/data_users/hkli/CFHT/w1234/original/source.list"

f = open(para_mode, "r")
contents = f.readlines()
f.close()

# initializing the directories "../astrometry, ../result, ../stamps"
cmd = "sh cfht_init.sh"
a = Popen(cmd, shell=True)
a.wait()

# run the three processes
process = 3
t = []
for i in range(process):
    t1 = time.clock()

    # change the PROCESS_stage
    contents[10] = "	parameter (PROCESS_stage=%d)\n"%i
    f = open(para_process, "w")
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
