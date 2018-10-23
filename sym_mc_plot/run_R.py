from subprocess import Popen

sources = ["pts", "ptsb", "dimmer", "dimmerm3"]

for source in sources:
    cmd = "mpirun -np 50 python R-factor.py %s"%source
    a = Popen(cmd, shell=True)
    a.wait()

