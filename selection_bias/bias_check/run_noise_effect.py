from subprocess import Popen

for i in range(4):
    cmd = "mpirun -np 30 ./noise_test %d"%i
    a = Popen(cmd, shell=True)
    a.wait()
