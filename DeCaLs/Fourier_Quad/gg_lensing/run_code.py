from sys import path,argv
path.append("/home/hklee/work/mylib")
from subprocess import Popen
import shutil
import hk_tool_box
import numpy


z = [[0.1, 2.0]]


dz = hk_tool_box.set_bin_log(0.001, 0.5, 31)
z_min, z_max = 0.15, 0.25

total_path = "/home/hklee/work/DECALS/DECALS_shear_catalog_v210729/gg_lensing"

# cmd = "mpirun -np 5 python prepare_cata.py prepare_foreground %.3f %.3f"%(z_min, z_max)
# a = Popen(cmd, shell=True)
# a.wait()


for iz in z:

    cmd = "mpirun -np 15 python prepare_cata.py prepare_background %.3f %.3f" % (iz[0], iz[1])
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "mpirun -np 1 python prepare_cata.py prepare_pdf"
    a = Popen(cmd, shell=True)
    a.wait()

    for tag, i_dz in enumerate(dz):

        cmd = "mpirun -np 40 ./GGL_cal %s cata/background cata/foreground 200 %.6f"%(total_path, i_dz)
        print(cmd)
        a = Popen(cmd, shell=True)
        a.wait()

        src = total_path + "/result/result.hdf5"
        dst = total_path + "/result/result_forez_%.2f_%.2f_dz_%d.hdf5"%(z_min, z_max, tag)
        shutil.move(src, dst)
