from subprocess import Popen
import time

data_type = ["noise_free", "noisy_cpp"]
for dt in data_type:
    for shear_tag in [2*i for i in range(15)]:
        t1 = time.time()
        cmd = "mpirun -np 30 python cal_chisq.py /mnt/perc/hklee/PDF_test/data %d %s"%(shear_tag, dt)
        print(cmd)
        a = Popen(cmd, shell=True)
        a.wait()
        t2 = time.time()
        print(dt,shear_tag, t2-t1)