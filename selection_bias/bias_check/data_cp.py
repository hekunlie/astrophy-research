from subprocess import Popen

total_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/selection_bias/bias_check/pts_sample/change_flux"

for i in range(5):
    src_path = "%s/imgs_%d/data"%(total_path, i)
    dst_path = "/mnt/ddnfs/data_users/hkli/bias_check/pow_noise_test/pts_sample/change_flux/imgs_%d"%i
    cmd = "scp -r %s hkli@202.120.32.231:%s"%(src_path, dst_path)
    print(cmd)
    a = Popen(cmd, shell=True)
    a.wait()
