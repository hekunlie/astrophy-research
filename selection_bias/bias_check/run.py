from subprocess import Popen

for i in range(17):
    cmd = "mpirun -np 40 python data_add.py /mnt/perc/hklee/bias_check/data_from_pi/new_data/data %d"%i
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "mpirun -np 40 ./calculate /mnt/perc/hklee/bias_check/data_from_pi/new_data/data/mix2 NF_add_CT_add_rpureCTest_rpsf_gg_%d mix 5"%i
    a = Popen(cmd, shell=True)
    a.wait()
