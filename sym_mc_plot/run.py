from subprocess import Popen


file_num = 1
# sources = ["dimmerm3", "dimmer", "pts", "ptsb"]
# resolution_factor = [0.7, 0.7, 0.33, 0.33]
sources = ["dimmerm3"]#, "pts", "ptsb"]
resolution_factor = [0.7]#, 0.33, 0.33]

cuts = ['sex_snr', 'mag_auto', "flux2"]
filter_name = ["sex2_2", "sex3_2","sex4_2"]
sig = ["2sig", "2sig", "2sig"]

# cuts = ['sex_snr', 'mag_auto', "flux2"]
# filter_name = ["sex2_1.5", "sex2_2", "sex3_1.5", "sex3_2", "sex4_1.5", "sex4_2"]
# sig = ["1.5sig", "2sig", "1.5sig", "2sig", "1.5sig", "2sig"]
for s, source in enumerate(sources):
    for f in range(len(filter_name)):
        for i in range(len(cuts)):
            if i < 2 or "sex2_" in filter_name[f]:
                cmd = "mpirun -np 14 python sym_mc_plot.py %s %d %s %s %s %.2f"%\
                      (cuts[i], file_num, filter_name[f], source, sig[f], resolution_factor[s])
                a = Popen(cmd, shell=True)
                a.wait()
