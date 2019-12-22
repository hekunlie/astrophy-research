from subprocess import Popen
import time


# total_path = "/mnt/perc/hklee/selection_bias/"
total_path = "/mnt/ddnfs/data_users/hkli/"

sources = ["pts_bright", "pts_dimmer"]

filter_name = ["sex2_1.5", "sex4_1.5", "sex2_2", "sex4_2", "sex2_4", "sex4_4"]
cuts = ["mag_true", "mag_auto", "rfactor", "flux2_ex1", "flux2_ex2", "flux2_ex3", "snr_sex"]


for s, source in enumerate(sources):
    for f in range(len(filter_name)):

        # for i in range(len(cuts)):
        #     t1 = time.time()
        #     cmd = "mpirun -np 50 /home/hkli/work/selection_bias/sym_mc_plot/cutoff %s %s %s"%\
        #           (total_path+source, filter_name[f], cuts[i])
        #     print("START: %s"%source)
        #     a = Popen(cmd, shell=True)
        #     a.wait()
        #     t2 = time.time()
        #     print("END: %s, %.2f sec"%(source, t2-t1))
        #
        #     print("Plot %s %s %s"%(source,filter_name[f],cuts[i]))
        #     cmd = " mpirun -np 10 python cutoff_plot.py %s%s %s %s"%(total_path, source,filter_name[f],cuts[i])
        #     a = Popen(cmd, shell=True)
        #     a.wait()

        print("Plot comparison")
        cmd = "python cutoff_plot_compare.py %s %s" % (total_path+source+"/result/cuts/sym", filter_name[f])
        a = Popen(cmd, shell=True)
        a.wait()