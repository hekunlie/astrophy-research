from subprocess import Popen
import time
import os
import shutil
from sys import argv
import logging

# specify the number of cores used
ncpus = 64

# specify it to where the para.inc is stored
para = "/home/hkli/work/fresh/src_new/para.inc"

# specify it to where the data are stored (don't forget the last "/")
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/w1234/original/"
source_list = data_path + "source.list"

# log
logger = logging.getLogger()
logger.setLevel(logging.INFO)

logfile = "/home/hkli/work/logs/cfht_log.dat"
lf = logging.FileHandler(logfile, 'w')
form = logging.Formatter('%(asctime)s - %(message)s')
lf.setFormatter(form)
logger.addHandler(lf)

# run the three processes
# one can specify the PROCESS_stage like [1,3], which comes from the commands inputs
# something likes "python fresh_run.py 1 2 3"
PROCESS_stage = []
nargv = len(argv)
if nargv == 1:
    print("No PROCESS_stage has been inputted.")
    exit()
else:
    for i in range(nargv):
        if 0 < i < 4:
            c = int(float(argv[i]))
            if c in [1, 2, 3]:
                PROCESS_stage.append(c)
            else:
                print("Wrong PROCESS_stage input: %s which should be one of [1, 2, 3]"%argv[i])
                exit()
    PROCESS_stage = sorted(PROCESS_stage)
print("The program will run the PROCESS_stage: ", PROCESS_stage)
time.sleep(3)
logger.info("PROCESS_stage: %s"%PROCESS_stage)

# read the para.inc
f = open(para, "r")
contents = f.readlines()
f.close()

files = os.listdir(data_path)
initial_target = ["/astrometry", "/stamps", "/result"]

t = [0, 0, 0]

for i in PROCESS_stage:
    t1 = time.time()
    logger.info("PROCESS_stage %d start"%i)
    # initialize the directories, "../astrometry" or " ../result" or "../stamps" (depending on the PROCESS_stage)
    # it will take a while...
    for file in files:
        if "w" in file:
            initial_path = data_path + file + initial_target[i-1]
            print("Deleting the %s"%initial_path)
            shutil.rmtree(initial_path, ignore_errors=True)
            os.mkdir(initial_path)
    logger.info("Deleted the %s's"%initial_target[i-1])

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
    logger.info("PROCESS_stage %d finish." %i)

print("PROCESS 1: %.2f H, PROCESS 2: %.2f H, PROCESS 3: %.2f H"%(t[0], t[1], t[2]))
