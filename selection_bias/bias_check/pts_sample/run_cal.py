from subprocess import Popen

# calculate shear
for ig in range(5):
    src_path = "/mnt/ddnfs/data_users/hkli/bias_check/pow_noise_test/pts_sample/change_flux/imgs_%d"%ig


    cmds = ["mpirun -np 10 %s/code/pdf_iter %s/data noisy_cpp 2000"%(src_path, src_path)]

    for cmd in cmds:
        a = Popen(cmd, shell=True)
        a.wait()
        print(cmd)

exit()
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
        "mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 cross_term noise_residual"%(src_path, src_path),
        "mpirun -np 10 %s/code/calculate_chisq %s/data 10 2000 noise_residual"%(src_path, src_path)]

for cmd in cmds:
    a = Popen(cmd, shell=True)
    a.wait()
    print(cmd)

