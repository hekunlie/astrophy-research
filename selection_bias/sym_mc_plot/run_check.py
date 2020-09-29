from subprocess import Popen

for i in range(20):
    cmd = "mpirun -np 6 python prepare_data.py /coma/hklee %d" %i
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "mpirun -np 15 python cutoff_cfht.py /coma/hklee/selection_bias %d" %i
    a = Popen(cmd, shell=True)
    a.wait()
