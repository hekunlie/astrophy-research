from subprocess import Popen
import os
import shutil
import time

with open("/home/hklee/work/catalog/Jesse_cata/hdf5/file_list","r") as f:
    files = f.readlines()

total_path = "/home/hklee/work/DECALS/DECALS_shear_catalog_v210729/gg_lensing"
sub_samples = [100, 200]

for fnm in files:
    if ".hdf5" in fnm:
        for j in range(len(sub_samples)):

            cmd = "mpirun -np 1 python prepare_cata.py prepare_foreground %s %d"%(fnm.split("\n")[0], sub_samples[j])
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)

            cmd = "mpirun -np 1 python prepare_cata.py prepare_pdf"
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)

            cmd = "mpirun -np 50 ./GGL_cal %s cata/background cata/foreground %d 0.2"%(total_path, sub_samples[j])
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)

            src = total_path + "/result/result.hdf5"
            dst = total_path + "/result_jesse/%s_%d_sub_sample.hdf5"%(fnm.split("_DECALS")[0], sub_samples[j])
            shutil.move(src, dst)
            #
            print("\nmove to %s\n"%dst)
            #src = total_path + "/result/result.hdf5"
            #dst = total_path + "/result/result_%d.hdf5"%j
            #shutil.move(src, dst)

            # print("\nmove to %s\n"%dst)
exit()
for j in range(3):
    for fnm in [files[1]]:
        if ".hdf5" in fnm:

            cmd = "mpirun -np 1 python prepare_cata_asym.py prepare_foreground %s"%fnm.split("\n")[0]
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)
            
            cmd = "mpirun -np 1 python prepare_cata_asym.py prepare_pdf"
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)

            cmd = "mpirun -np 50 ./GGL_cal /home/hklee/work/DECALS/gg_lensing 100 0.2"
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)

            src = total_path + "/result/result.hdf5"
            dst = total_path + "/result_jesse/guess_test/asym/%s_%d.hdf5"%(fnm.split("_DECALS")[0],j)
            shutil.move(src, dst)
