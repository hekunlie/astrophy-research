from subprocess import Popen
from sys import argv

cmd = argv[1]

if cmd == "simu":
    for i in range(1, 6):
        cmd = "mpirun -np 50 python /mnt/ddnfs/data_users/hkli/bias_check/test/check_ct_property.py simu %d 46212"%i
        a = Popen(cmd, shell=True)
        a.wait()
        print(cmd)
else:
    cmd = "mpirun -np 5 python /mnt/ddnfs/data_users/hkli/bias_check/test/check_ct_property.py chisq"
    a = Popen(cmd, shell=True)
    a.wait()
    print(cmd)