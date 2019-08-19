import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
import h5py
from mpi4py import MPI
import numpy
import time
from subprocess import Popen
import warnings
warnings.filterwarnings('error')



# The new Fourier_Quad catalog differs from the old version!!!
# collect: collect the data from the files of each field. It creates the "fourier_cata.hdf5" in
#           the parent directory of the one contain the field catalog.
#           If the catalog file doesn't exist, run it firstly !!!.
#           It will add the redshift parameters from CFHT catalog into the finial catalog.


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

data_path = "/mnt/perc/hklee/CFHT/catalog/fourier_cata_new/"
raw_cata_path = data_path + "raw_cata_new/"

dicts, fields = tool_box.field_dict(data_path + "nname.dat")

my_field = tool_box.allot(fields, cpus)[rank]

chip_num = 36

for field_nm in my_field:
    field_path = raw_cata_path + "%s/"%field_nm
    files = os.listdir(field_path)

    chip_exps = []
    for nm in files:
        if ".dat" in nm:
            exp_nm = nm.split("p")[0]
            if exp_nm not in chip_exps:
                chip_exps.append(exp_nm)
    chip_exps.sort()

    file_count = 0
    for exp_nm in chip_exps:
        for i in range(1,chip_num+1):
            chip_nm = "%sp_%d_shear.dat"%(exp_nm, i)
            chip_path = field_path + chip_nm
            if os.path.exists(chip_path):
                try:
                    temp = numpy.loadtxt(chip_path, skiprows=1)
                    if file_count == 0:
                        data = temp
                    else:
                        data = numpy.row_stack((data, temp))
                    file_count += 1
                except:
                    file_size = os.path.getsize(chip_path)/1024.
                    print("Empty: %s (%.3f KB)"%(chip_nm, file_size))
            else:
                print("Can't find %d"%chip_nm)
    if file_count > 0:
        final_path = data_path + "%s/%s_shear_raw.cat"%(field_nm, field_nm)
        numpy.savetxt(final_path, data)
        h5f = h5py.File(final_path,"w")
        h5f["/data"] = data
        h5f.close()
