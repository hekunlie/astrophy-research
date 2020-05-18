from subprocess import Popen
from sys import argv

total_path = ["/mnt/perc/hklee/bias_check/data_from_pi/pow_noise_test/noise_free",
                "/mnt/perc/hklee/bias_check/data_from_pi/pow_noise_test/noisy"]

for src_path in total_path:
    for dg in ["0.04", "0.01","0.005","0.002"]:
        cmd = "mpirun -np 10 ./new_pdf " + src_path + " data_%d_epsf.hdf5 4 100 " + dg
        print(cmd)
        a = Popen(cmd, shell=True)
        a.wait()