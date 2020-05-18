import os
from sys import path, argv
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
# path.append('%s/work/mylib/' % my_home)
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
import numpy
import h5py
import tool_box


result = numpy.zeros((2,4))
for j in range(1):
    data_path = "E:/data/new_pdf/pts_sample/imgs_%d"%j
    h5f = h5py.File(data_path + "/shear.hdf5", "r")
    g1 = h5f["/g1"][()]
    g2 = h5f["/g2"][()]
    h5f.close()
    h5f = h5py.File(data_path + "/shear_result_pdf_iter_noisy_cpp.hdf5", "r+")
    mg1 = h5f["/new_PDF/g1"][()]
    mg2 = h5f["/new_PDF/g2"][()]
    h5f.close()
    print(mg1.shape)
    for i in range(3):
        mc1 = tool_box.data_fit(g1, mg1[:,i*6],mg1[:,i*6+1])
        mc2 = tool_box.data_fit(g2, mg2[:,i*6],mg2[:,i*6+1])
        result[0] = mc1[0]-1, mc1[1], mc1[2], mc1[3]
        result[1] = mc2[0]-1, mc2[1], mc2[2], mc2[3]

    # h5f["/new_PDF_sym_mc"] = result


        print(result)




