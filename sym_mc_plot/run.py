from subprocess import Popen

# sources = ["dimmerm3", "dimmer", "pts", "ptsb"]
# resolution_factor = [0.7, 0.7, 0.33, 0.33]
sources = ["dimmerm3"]#, "pts", "ptsb"]
resolution_factor = [0.7]#, 0.33, 0.33]

cuts = ["snr_auto", 'sex_snr', 'mag_auto', "flux2","flux_alt"]
filter_name = ["sex2_1.5", "sex3_1.5","sex4_1.5","sex2_2", "sex3_2","sex4_2"]
sig = ["1.5sig","1.5sig","1.5sig","2sig", "2sig", "2sig"]

# cuts = ['sex_snr', 'mag_auto', "flux2"]
# filter_name = ["sex2_1.5", "sex2_2", "sex3_1.5", "sex3_2", "sex4_1.5", "sex4_2"]
# sig = ["1.5sig", "2sig", "1.5sig", "2sig", "1.5sig", "2sig"]
for s, source in enumerate(sources):
    for f in range(len(filter_name)):
        for i in range(len(cuts)):
            if i < 3 or "sex2_" in filter_name[f]:
                cmd = "mpirun -np 10 python sym_mc_plot.py %s %s %s %s %.2f"%\
                      (cuts[i],  filter_name[f], source, sig[f], resolution_factor[s])
                a = Popen(cmd, shell=True)
                a.wait()
