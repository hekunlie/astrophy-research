from subprocess import Popen
import time
from sys import argv
import os

sources = [argv[1]]
resolution_factor = [0.7]
while True:
    if os.path.exists("finish.dat"):
        break
if argv[2] == "sex":
    filter_name = ["sex2_1.5", "sex2_2","sex3_1.5","sex3_2", "sex4_1.5","sex4_2"]
    cuts = ["mag_auto", "snr_auto", 'sex_snr']
elif argv[2] == "f1":
    filter_name = ["sex2_1.5", "sex3_1.5","sex4_1.5"]
    cuts = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]
elif argv[2] == "f2":
    filter_name = ["sex2_2", "sex3_2","sex4_2"]
    cuts = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]
else:
    print("Wrong!!! %s"%argv[2])
    exit()
for s, source in enumerate(sources):
    for f in range(len(filter_name)):
        for i in range(len(cuts)):
            t1 = time.time()
            sig = float(filter_name[f].split("_")[-1])
            cmd = "mpirun -np 10 python sym_mc_plot.py %s %s %s %.1fsig %.2f"%\
                  (cuts[i],  filter_name[f], source, sig, resolution_factor[s])
            print("START:",cmd.split("py")[-1])
            a = Popen(cmd, shell=True)
            a.wait()
            t2 = time.time()
            print("END:", cmd.split("py")[-1], t2-t1)