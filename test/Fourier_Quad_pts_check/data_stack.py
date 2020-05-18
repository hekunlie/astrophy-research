import numpy
import h5py
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

total_path = "/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/code_test/new_pdf_test/data/epsf"

data = numpy.zeros((7000000, 8))
sub_row = 1400000

for tag, kk in enumerate([1.0, 1.5, 2.0, 2.5, 3.0]):
    h5f = h5py.File(total_path + "/result_%.1f/data_%d_noise_free.hdf5"%(kk,rank), "r")
    data[tag*sub_row:(tag+1)*sub_row] = h5f["/data"][()][:sub_row]
    h5f.close()

h5f = h5py.File(total_path+"/result_stack/data_%d_noise_free.hdf5"%rank, "w")
h5f["/data"] = data
h5f.close()

for tag, kk in enumerate([1.0, 1.5, 2.0, 2.5, 3.0]):
    h5f = h5py.File(total_path + "/result_%.1f/data_%d_noisy.hdf5"%(kk,rank), "r")
    data[tag*sub_row:(tag+1)*sub_row] = h5f["/data"][()][:sub_row]
    h5f.close()

h5f = h5py.File(total_path+"/result_stack/data_%d_noisy.hdf5"%rank, "w")
h5f["/data"] = data
h5f.close()