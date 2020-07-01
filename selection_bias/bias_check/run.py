from subprocess import Popen

parent_path = "/mnt/perc/hklee/bias_check/data_from_pi/new_data/data"
for j in range(4):
    for i in range(8):
        cmd = "mpirun -np 20 python data_add.py %s %d %d"%(parent_path, i, j)
        a = Popen(cmd, shell=True)
        a.wait()

        cmd = "mpirun -np 40 ./calculate %s/mix NF_add_CT_add_rpureCTest_rpsf_%d_%d mix 5"%(parent_path,j,i)
        a = Popen(cmd, shell=True)
        a.wait()

        print(j,i)
