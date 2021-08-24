from subprocess import Popen
import shutil


z = [[0.25,0.35], [0.35,0.45], [0.45,0.55], [0.55,0.65], [0.65,0.75], [0.75,0.85], [0.85,0.95],[0.95, 1.05], [1.05, 1.15]]



total_path = "/home/hklee/work/DECALS/DECALS_shear_catalog_old/gg_lensing"

cmd = "mpirun -np 5 python prepare_cata.py prepare_foreground"
a = Popen(cmd, shell=True)
a.wait()

for iz in z:
    cmd = "mpirun -np 15 python prepare_cata.py prepare_background %.2f %.2f"%(iz[0], iz[1])
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "mpirun -np 1 python prepare_cata.py prepare_pdf"
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "mpirun -np 50 ./GGL_cal %s cata/background cata/foreground 200 0.0"%total_path
    a = Popen(cmd, shell=True)
    a.wait()

    src = total_path + "/result/result.hdf5"
    dst = total_path + "/result/result_%.2f_%.2f_dz_0.1.hdf5"%(iz[0], iz[1])
    shutil.move(src, dst)