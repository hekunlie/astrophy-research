import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
import numpy
from sys import path, argv
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
path.append('%s/work/mylib/' % my_home)
import h5py
import tool_box


file_path = argv[1]

h5f = h5py.File(file_path,"r")
data = h5f["/field"][()]
expo_label = data[:,-1]
print(expo_label)
min_label, max_label = expo_label.min(), expo_label.max()
expo_num = int(max_label-min_label+1)
expo_num_ = h5f["/expos_num"][()][0]
expo_label = (expo_label - min_label).astype(dtype=numpy.intc)
print(expo_num,expo_num_)
print(expo_label)

for i in range(expo_num):
    idx = expo_label == i

    expo_name = h5f["/expo_%d"%i].attrs["exposure_name"]
    print(expo_name)
    expo_data_h5 = h5f["/expo_%d"%i][()][:,:5]
    diff = data[idx][:,:5] - expo_data_h5
    print(diff.min(), diff.max())

h5f.close()