import h5py
from sys import argv

# delete the specific data set in a hdf5 file
data_path = argv[1]
set_name = argv[2]

h5f = h5py.File(data_path)
del h5f[set_name]
h5f.close()
print("Delete: %s"%set_name)
