import h5py
import numpy
from sys import argv

# choose a subset of foreground catalog
area_id = int(argv[1])
num = int(argv[2])

h5f = h5py.File("./data/w_%d.hdf5"%area_id,"r")
nms = list(h5f.keys())
datas = []
for nm in nms:
    datas.append(h5f[nm].value)
h5f.close()

total_num = datas[0].shape[0]
ch = numpy.random.choice(numpy.arange(total_num),num, False)

h5f_new = h5py.File("./data/w_%d_sub.hdf5"%area_id,"w")
for i in range(len(nms)):
    h5f_new[nms[i]] = datas[i][ch]
h5f_new.close()
print("Chose %d galaxies"%num)