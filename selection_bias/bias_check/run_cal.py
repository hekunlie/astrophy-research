from subprocess import Popen


data_row = [5000, 1000, 1000, 1000, 1000, 1000]
for i in range(6,12):
    cmd = "mpirun -np 20 ./calculate /mnt/perc/hklee/bias_check/result/data_collection/%d %d 8"%(i, data_row[i-6])
    a = Popen(cmd, shell=True)
    a.wait()
