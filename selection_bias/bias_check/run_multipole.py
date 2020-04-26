from subprocess import Popen

src_nm = ["galsim_sample", "pts_sample"]
data_path = "/mnt/ddnfs/data_users/hkli/bias_check/pow_noise_test"
for src in src_nm:
    for i in range(5):
        cmd = "mpirun -np 10 python multi_pole_fit.py %s/%s/change_flux/imgs_%d/data 200 80"%(data_path, src, i)
        a = Popen(cmd, shell=True)
        a.wait()
