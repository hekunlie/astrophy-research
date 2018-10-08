from subprocess import Popen


file_num = 5
cuts = ['sex_snr', 'mag_auto']
filter_name = ["sex2_1.5", "sex2_2", "sex3_1.5", "sex3_2", "sex4_1.5", "sex4_2"]

for f in range(len(filter_name)):
    for i in range(len(cuts)):
        cmd = "mpirun -np 14 python sym_mc_plot.py %s %d %s"%(cuts[i], file_num, filter_name[f])
        a = Popen(cmd, shell=True)
        a.wait()