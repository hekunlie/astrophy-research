from subprocess import Popen

for i in range(1, 7):
    for j in range(18):
        cmd = "python show_mass_map.py %d 1.0 %d"%(j,i)
        a = Popen(cmd, shell=True)
        a.wait()
