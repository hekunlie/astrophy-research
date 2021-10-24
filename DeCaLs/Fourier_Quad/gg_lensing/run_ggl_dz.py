from subprocess import Popen
import os
import shutil
import time
import h5py
import numpy

total_path = "/home/hklee/work/DECALS/DECALS_shear_catalog_v210729/gg_lensing"

rand_fore_path = "/home/hklee/work/catalog/Yang_group/DESI_CLUSTER_DR9"

# cmd = "mpirun -np 30 python prepare_cata.py prepare_background 0 3"
# a = Popen(cmd, shell=True)
# a.wait()
# time.sleep(2)

dz = 0.005

fore_z = [0.2 + i*0.1 for i in range(8)]

for iz in fore_z:

    cmd = "mpirun -np 5 python prepare_cata.py prepare_foreground %.4f %.4f DESIDR9_NGC_SGC_group_DECALS_overlap.hdf5"%(iz, iz+dz)
    a = Popen(cmd, shell=True)
    a.wait()
    time.sleep(2)

    cmd = "mpirun -np 55 ./ggl_dz_hist %s cata/background cata/foreground 200 result_len_%.4f"%(total_path, iz)
    a = Popen(cmd, shell=True)
    a.wait()
    time.sleep(2)

    src = total_path + "/result/result_len_%.4f.hdf5"%iz
    dst = total_path + "/result_dz/result_len_%.4f.hdf5"%iz
    shutil.move(src, dst)
    #
    print("\nmove to %s\n" % dst)


    #  prepare the random catalog
    h5f = h5py.File(rand_fore_path + "/rand_fore_candidate.hdf5", "r")
    data = h5f["/data"][()]
    h5f.close()
    rand_num = data.shape[0]

    h5f = h5py.File(rand_fore_path + "/rand_fore.hdf5", "w")
    data[:, 3] = numpy.random.uniform(iz, iz+dz, rand_num).astype(numpy.float32)
    data[:,4] = 14.1
    h5f["/data"] = data
    h5f.close()

    cmd = "mpirun -np 5 python prepare_cata.py prepare_foreground %.4f %.4f rand_fore.hdf5"%(iz, iz+dz)
    a = Popen(cmd, shell=True)
    a.wait()
    time.sleep(2)

    cmd = "mpirun -np 55 ./ggl_dz_hist %s cata/background cata/foreground 200 result_rand_%.4f" % (total_path, iz)
    a = Popen(cmd, shell=True)
    a.wait()
    time.sleep(2)

    src = total_path + "/result/result_rand_%.4f.hdf5"%iz
    dst = total_path + "/result_dz/result_rand_%.4f.hdf5"%iz
    shutil.move(src, dst)
    #
    print("\nmove to %s\n" % dst)