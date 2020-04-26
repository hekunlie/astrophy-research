from subprocess import Popen
from sys import argv

parent_path = "/mnt/ddnfs/data_users/hkli/bias_check/pow_noise_test/change_flux"
code_path = parent_path +"/imgs_5/code/simu"

seeds = [1243,346547,215632,78996,123566,986943,456012]
for i in range(5,6):
    data_path = "%s/imgs_%d"%(parent_path, i)
    cmd = "mpirun -np 100 %s %s %d %d 2000 44"%(code_path, data_path, seeds[i], i)
    a = Popen(cmd, shell=True)
    a.wait()


# calculate shear
src_path = "/mnt/ddnfs/data_users/hkli/bias_check/pow_noise_test/change_flux/imgs_5"

cmds = ["mpirun -np 10 %s/code/calculate_s %s/data noise_free 10 2000"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_s %s/data noisy_cpp 10 2000"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_m %s/data N_add_rCT noisy_cpp cross_term 1 10 2000"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_m %s/data N_add_rCT_est noisy_cpp cross_term_est 1 10 2000"%(src_path, src_path)]

for cmd in cmds:
    a = Popen(cmd, shell=True)
    a.wait()
    print(cmd)

# calculate chisq
cmds = ["mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 noise_free"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 noisy_cpp"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 cross_term"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 cross_term_sqrt"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 cross_term noise_residual"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 noise_residual"%(src_path, src_path)]

for cmd in cmds:
    a = Popen(cmd, shell=True)
    a.wait()
    print(cmd)
