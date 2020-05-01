from subprocess import Popen
from sys import argv

src_nm = [argv[1]]
shear_scale = argv[2]
data_path = "/mnt/ddnfs/data_users/hkli/bias_check/pow_noise_test"
for src in src_nm:
    for i in range(5):
        cmd = "mpirun -np 10 python multi_pole_fit.py %s/%s/change_flux/e_pdf/imgs_%d/data 200 60 %s"%(data_path, src, i,shear_scale)
        a = Popen(cmd, shell=True)
        a.wait()
