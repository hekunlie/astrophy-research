from subprocess import Popen
import os
import shutil
import time


para_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_10"
code_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/code"
for isig in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]:

    cmd = "mpirun -np 55 python %s/NFW_shear_map_continue_source_plane.py 0 %.2f 4000000 3"%(para_path, isig)
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "python %s/stack_file.py %s 4000000"%(code_path, para_path)
    a = Popen(cmd, shell=True)
    a.wait()

    src = "%s/params/segment_sheared_para.hdf5"%para_path
    dst = "%s/params/segment_sheared_para_lenz_%.2f.hdf5"%(para_path, isig)
    os.rename(src, dst)
exit()

para_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_9"
code_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/code"
for isig in [0.05, 0.1, 0.15]:

    cmd = "mpirun -np 55 python %s/NFW_shear_map_continue_source_plane.py 0 2000000 3 %.2f"%(para_path, isig)
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "python %s/stack_file.py %s 2000000"%(code_path, para_path)
    a = Popen(cmd, shell=True)
    a.wait()

    cmd = "python %s/segment_file.py "%code_path
    a = Popen(cmd, shell=True)
    a.wait()

    src = "%s/params/segment_sheared_para.hdf5"%para_path
    dst = "%s/params/segment_sheared_para_decals_Pz_errsig_z%.2f.hdf5"%(para_path, isig)
    os.rename(src, dst)

cmd = "mpirun -np 3 python %s/zerr_test.py"%code_path
a = Popen(cmd, shell=True)
a.wait()
