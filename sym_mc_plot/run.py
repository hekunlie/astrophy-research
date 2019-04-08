from subprocess import Popen
import time
from sys import argv
import os

sources = [argv[1]]
# cut_criterion = argv[2]
# resolution_factor = [0.7]
# while True:
#     if os.path.exists("finish.dat"):
#         break
# if cut_criterion == "s1":
#     filter_name = ["sex2_1.5", "sex3_1.5", "sex4_1.5"]
#     cuts = ["mag_auto", "snr_auto", 'sex_snr']
#
# elif cut_criterion == "s2":
#     filter_name = ["sex2_2", "sex3_2", "sex4_2"]
#     cuts = ["mag_auto", "snr_auto", 'sex_snr']
#
# elif cut_criterion == "s3":
#     filter_name = ["sex2_4", "sex3_4", "sex4_4"]
#     cuts = ["mag_auto", "snr_auto", 'sex_snr']
#
# elif cut_criterion == "f1":
#     filter_name = ["sex2_1.5", "sex3_1.5", "sex4_1.5"]
#     cuts = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]
#
# elif cut_criterion == "f2":
#     filter_name = ["sex2_2", "sex3_2", "sex4_2"]
#     cuts = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]
#
# elif cut_criterion == "f3":
#     filter_name = ["sex2_4", "sex3_4", "sex4_4"]
#     cuts = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]
#
# else:
#     print("Wrong!!! %s"%argv[2])
#     exit()

filter_name = ["sex2_1.5", "sex3_1.5", "sex4_1.5","sex2_2", "sex3_2", "sex4_2","sex2_4", "sex3_4", "sex4_4"]
cuts = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]


for s, source in enumerate(sources):
    for f in range(len(filter_name)):
        for i in range(len(cuts)):
            t1 = time.time()
            sig = float(filter_name[f].split("_")[-1])
            # cmd = "mpirun -np 10 python sym_mc_plot.py %s %s %s %.1fsig %.2f"%\
            #       (cuts[i],  filter_name[f], source, sig, resolution_factor[s])
            cmd = "mpirun -np 10 /home/hkli/work/cpp/cutoff %s %s %.1f %s"%\
                  (source, filter_name[f], sig, cuts[i])
            print("START:",cmd.split("py")[-1])
            a = Popen(cmd, shell=True)
            a.wait()
            t2 = time.time()
            print("END:", cmd.split("py")[-1], t2-t1)