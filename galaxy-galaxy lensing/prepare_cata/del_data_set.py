import h5py
from sys import argv

# delete the specific data set in a hdf5 file
data_path = argv[1]
num = len(argv)

h5f = h5py.File(data_path)
for i in range(2,num):
    del h5f[argv[i]]
    print("Delete: %s" % argv[i])
h5f.close()

