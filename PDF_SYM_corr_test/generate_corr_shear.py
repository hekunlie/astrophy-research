import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
import numpy
import tool_box
import time


cor_gg_len = int(argv[1])
cmd = argv[2]

chi_guess_bin = [1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9]
chi_guess_num = len(chi_guess_bin)

shear_test = numpy.linspace(-0.04,0.04, chi_guess_num)
mean = [0, 0]

for i in range(chi_guess_num):

    if cmd == "corr":
        cov = [[numpy.abs(chi_guess_bin[i] * 2), chi_guess_bin[i]],
               [chi_guess_bin[i], numpy.abs(chi_guess_bin[i] * 2)]]

        gg1 = tool_box.rand_gauss2n(cor_gg_len, mean, cov).astype(dtype=numpy.float32)
        gg2 = tool_box.rand_gauss2n(cor_gg_len, mean, cov).astype(dtype=numpy.float32)
    else:
        gg1 = numpy.ones((2,cor_gg_len))*shear_test[i]
        gg2 = numpy.ones((2,cor_gg_len))*shear_test[i]*(-1)

    for j in range(2):
        print(i,j,2*i+j)
        h5f = h5py.File("/mnt/perc/hklee/PDF_test/data/shear_%d.hdf5"%int(2*i+j),"w")
        h5f["/g1"] = gg1[j]
        h5f["/g2"] = gg2[j]
        h5f.close()
