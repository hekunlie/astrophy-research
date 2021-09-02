from subprocess import Popen
import shutil


z = [[0.1, 1.5]]

dz = [ 0.2, 0.25, 0.3, 0.35]#0.05, 0.1, 0.15,
z_min, z_max = 0.15, 0.25

total_path = "/home/hklee/work/DECALS/DECALS_shear_catalog_v210729/gg_lensing"

cmd = "mpirun -np 5 python prepare_cata.py prepare_foreground %.2f %.2f"%(z_min, z_max)
a = Popen(cmd, shell=True)
a.wait()


for iz in z:

    cmd = "mpirun -np 15 python prepare_cata.py prepare_background %.2f %.2f" % (iz[0], iz[1])
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "mpirun -np 1 python prepare_cata.py prepare_pdf"
    a = Popen(cmd, shell=True)
    a.wait()

    for i_dz in dz:

        cmd = "mpirun -np 50 ./GGL_cal %s cata/background cata/foreground 200 %.2f"%(total_path, i_dz)
        a = Popen(cmd, shell=True)
        a.wait()

        src = total_path + "/result/result.hdf5"
        dst = total_path + "/result/result_forez_%.2f_%.2f_dz_%.2f.hdf5"%(z_min, z_max, i_dz)
        shutil.move(src, dst)
