from subprocess import Popen
import os
import shutil
import time

cat_path = "/home/hklee/work/catalog/Jesse_cata/hdf5"
cats = os.listdir(cat_path)

total_path = "/home/hklee/work/DECALS/gg_lensing"
for fnm in cats:
    if ".hdf5" in fnm:

        cmd = "mpirun -np 1 python prepare_cata.py prepare_foreground %s"%fnm
        a = Popen(cmd, shell=True)
        a.wait()
        time.sleep(2)

        cmd = "mpirun -np 1 python prepare_cata.py prepare_pdf"
        a = Popen(cmd, shell=True)
        a.wait()
        time.sleep(2)

        cmd = "mpirun -np 20 ./GGL_cal /home/hklee/work/DECALS/gg_lensing 100 0.2"
        a = Popen(cmd, shell=True)
        a.wait()
        time.sleep(2)

        src = total_path + "/result/result.hdf5"
        dst = total_path + "/result_jesse/%s.hdf5"%fnm.split("_DECALS")[0]
        shutil.move(src, dst)
