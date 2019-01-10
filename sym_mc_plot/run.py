from subprocess import Popen


sources = ["debug"]#, "pts", "ptsb"]
resolution_factor = [0.7]#, 0.33, 0.33]

# filter_name = ["sex2_1.5", "sex2_2","sex3_1.5","sex3_2", "sex4_1.5","sex4_2"]
# cuts = ["flux2", "flux_alt", "mag_auto", "snr_auto", 'sex_snr']

filter_name = ["sex2_1.5", "sex2_2"]
cuts = ["flux2", "flux_alt", "flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4","flux2_ex5"]
#cuts = ["flux2_ex2", "flux2_ex3", "flux2_ex4"]
#cuts = ["flux2_ex5", "flux_ex1"]#, "flux2_ex4"]
# cuts = ["flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]
# cuts = ["flux_ex1", "flux_ex2", "flux_ex3", "flux_ex4", "flux_ex5"]
#cuts = ["snr_ex1", "snr_ex2", "snr_ex3", "snr_ex4", "snr_ex5"]
for s, source in enumerate(sources):
    for f in range(len(filter_name)):
        for i in range(len(cuts)):
            # if i > 1 or "sex2_" in filter_name[f]:
            if "sex2_" in filter_name[f]:
                sig = float(filter_name[f].split("_")[-1])
                cmd = "mpirun -np 10 python sym_mc_plot.py %s %s %s %.1fsig %.2f"%\
                      (cuts[i],  filter_name[f], source, sig, resolution_factor[s])
                print(cmd.split("py")[-1])
                a = Popen(cmd, shell=True)
                a.wait()
