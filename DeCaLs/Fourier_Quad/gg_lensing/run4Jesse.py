from subprocess import Popen
import os
import shutil
import time

with open("/home/hklee/work/catalog/Jesse_cata/hdf5/file_list","r") as f:
    files = f.readlines()

total_path = "/home/hklee/work/DECALS/gg_lensing"
for j in range(5):
    for fnm in [files[0]]:
        if ".hdf5" in fnm:

            # cmd = "mpirun -np 1 python prepare_cata.py prepare_foreground %s"%fnm.split("\n")[0]
            # a = Popen(cmd, shell=True)
            # a.wait()
            # time.sleep(2)
            #
            cmd = "mpirun -np 1 python prepare_cata.py prepare_pdf"
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)

            cmd = "mpirun -np 50 ./GGL_cal /home/hklee/work/DECALS/gg_lensing 50 0.2"
            a = Popen(cmd, shell=True)
            a.wait()
            time.sleep(2)

            # src = total_path + "/result/result.hdf5"
            # dst = total_path + "/result_jesse/%s.hdf5"%fnm.split("_DECALS")[0]
            # shutil.move(src, dst)
            #
            # print("\nmove to %s\n"%dst)
            src = total_path + "/result/result.hdf5"
            dst = total_path + "/result/result_%d.hdf5"%j
            shutil.move(src, dst)

            print("\nmove to %s\n"%dst)