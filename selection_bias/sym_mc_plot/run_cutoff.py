from subprocess import Popen
import time
from sys import argv
import os

# sources = [argv[1]]
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

# total_path = "/mnt/perc/hklee/selection_bias/"
total_path = "/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/"

sources = ["galsim_bright", "galsim_dimmer", "pts_bright", "pts_dimmer"]
# sources = ["pts_dimmer"]
filter_name = ["sex2_1.5", "sex4_1.5", "sex2_2", "sex4_2", "sex2_4", "sex4_4"]
cuts = ["mag_auto", "rfactor", "flux2_ex1", "flux2_ex2", "flux2_ex3", "snr_sex","snr_auto"]


for s, source in enumerate(sources):
    for f in range(len(filter_name)):
        for i in range(len(cuts)):
            t1 = time.time()
            # cmd = "mpirun -np 30 /home/hklee/work/selection_bias/sym_mc_plot/cutoff %s %s %s"%\
            #       (total_path+source, filter_name[f], cuts[i])
            cmd = "mpirun -np 50 /home/hkli/work/selection_bias/sym_mc_plot/cutoff %s %s %s"%\
                  (total_path+source, filter_name[f], cuts[i])
            print("START: %s"%source)
            a = Popen(cmd, shell=True)
            a.wait()
            t2 = time.time()
            print("END: %s, %.2f sec"%(source, t2-t1))