from subprocess import Popen


for i in range(1,6):
    cmd = "mpirun -np 20 ./calculate /mnt/perc/hklee/bias_check/result/datas/%d 1000 8"%i
    a = Popen(cmd, shell=True)
    a.wait()
