from subprocess import Popen
import os

total_src_path = "/lustre/home/acct-phyzj/phyzj/CFHT/i"
# total_dst_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/fourier_cata"
total_dst_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/CFHT"

with open("nname_.dat", "r") as f:
    fields = f.readlines()

for nm in fields:
    if "w" in nm:
        try:
            nm = nm.split("\n")[0]
            dst_path = total_dst_path + "/%s"%nm
            if not os.path.exists(dst_path):
                os.makedirs(dst_path)

            src_path = "%s/%s/result"%(total_src_path, nm)

            # cmd = "sc -r %s hkli@202.120.32.231:%s"%(src_path, dst_path)
            cmd = "cp -r %s %s/"%(src_path, dst_path)
            print(cmd)
            a = Popen(cmd, shell=True)
            a.wait()
        except:
            pass
