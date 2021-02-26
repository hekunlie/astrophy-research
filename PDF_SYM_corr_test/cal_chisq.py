import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append("%s/work/mylib/"% my_home)
import h5py
import numpy
from mpi4py import MPI
# import tool_box
# from plot_tool import Image_Plot
import tool_box
from Fourier_Quad import Fourier_Quad
import time


def get_chisq(count, bin_num):
    ny, nx = int(bin_num / 2), int(bin_num / 2)
    arr_1 = count[0:ny, 0:nx][:, range(ny - 1, -1, -1)]
    arr_2 = count[0:ny, nx:2 * nx]
    arr_3 = count[ny:2 * ny, 0:nx][range(ny - 1, -1, -1)][:, range(nx - 1, -1, -1)]
    arr_4 = count[ny:2 * ny, nx:2 * nx][range(ny - 1, -1, -1)]
    return 0.5 * numpy.sum(((arr_2 + arr_3 - arr_1 - arr_4) ** 2) / (arr_1 + arr_2 + arr_3 + arr_4))



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

t1 = time.time()

data_path = argv[1]
tag_1 = int(argv[2])
tag_2 = tag_1 + 1
data_type = argv[3]

# chi guess bin for PDF_SYM
# chi_guess_num = 350
# inv = [chi_guess_num-1-i for i in range(chi_guess_num)]
# chi_guess_bin_p = tool_box.set_bin_log(3.*10**(-10), 6.*10**(-4), chi_guess_num).astype(numpy.float32)
# chi_guess_bin = numpy.zeros((2*chi_guess_num, ), dtype=numpy.float32)
# chi_guess_bin[:chi_guess_num] = -chi_guess_bin_p[inv]
# chi_guess_bin[chi_guess_num:] = chi_guess_bin_p
# chi_guess_bin = numpy.sort(chi_guess_bin)
#
# chi_guess_num = int(chi_guess_num*2)

chi_guess_num = 100
inv = [chi_guess_num-1-i for i in range(chi_guess_num)]
xi_cache = numpy.load(data_path+"/xi_cache_%s.npz"%data_type)["arr_0"]
signal_p = xi_cache[0,int(tag_1/2)]
signal_p_err = xi_cache[1,int(tag_1/2)]*10
chi_guess_bin = numpy.linspace(signal_p-signal_p_err, signal_p + signal_p_err, chi_guess_num).astype(numpy.float32)


mg_bin_num = 10


total_count_1 = numpy.zeros((chi_guess_num, mg_bin_num, mg_bin_num))
total_count_2 = numpy.zeros((chi_guess_num, mg_bin_num, mg_bin_num))


if rank == 0:
    print(tag_1, tag_2, data_path)

itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = 2*chi_guess_num*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
chisq_arr = numpy.ndarray(buffer=buf1, dtype='d', shape=(2,chi_guess_num))


task_list = [i for i in range(chi_guess_num)]
my_task_list = tool_box.alloc(task_list, cpus)[rank]


if rank == 0:
    if not os.path.exists(data_path + "/mg_bin_%s.hdf5"%data_type):
        fq = Fourier_Quad(12, 12412)
        h5f = h5py.File(data_path + "/data_%d_%s.hdf5"%(tag_1,data_type),"r")
        shear_data = h5f["/data_0"][()]
        h5f.close()
        mg_bin = fq.set_bin(shear_data[:,0], mg_bin_num)
        h5f = h5py.File(data_path + "/mg_bin_%s.hdf5"%data_type,"w")
        h5f["/data"] = mg_bin
        h5f.close()
        print(mg_bin)
comm.Barrier()

h5f = h5py.File(data_path + "/mg_bin_%s.hdf5" % data_type, "r")
mg_bin = h5f["/data"][()]
h5f.close()
bins = [mg_bin, mg_bin]

mean = [0, 0]

h5f_1 = h5py.File(data_path + "/data_%d_%s.hdf5"%(tag_1, data_type),"r")

h5f_2 = h5py.File(data_path + "/data_%d_%s.hdf5"%(tag_2, data_type),"r")

section_num = len(list(h5f_1.keys()))

for tag in range(section_num):

    shear_data_1 = h5f_1["/data_%d" % tag][()]
    shear_data_2 = h5f_2["/data_%d" % tag][()]
    cor_gg_len = shear_data_1.shape[0]

    for i in my_task_list:

        cov = [[numpy.abs(chi_guess_bin[i] * 2), -chi_guess_bin[i]],
               [-chi_guess_bin[i], numpy.abs(chi_guess_bin[i] * 2)]]

        gg1 = tool_box.rand_gauss2n(cor_gg_len, mean, cov).astype(dtype=numpy.float32)

        mg1_1 = shear_data_1[:,0] - gg1[0]*(shear_data_1[:,2] + shear_data_1[:,3])
        mg1_2 = shear_data_2[:,0] - gg1[1]*(shear_data_2[:,2] + shear_data_2[:,3])

        count_1 = numpy.histogram2d(mg1_1, mg1_2, bins)[0]
        total_count_1[i] = total_count_1[i] + count_1


        mg2_1 = shear_data_1[:,1] - gg1[0]*(shear_data_1[:,2] - shear_data_1[:,3])
        mg2_2 = shear_data_2[:,1] - gg1[1]*(shear_data_2[:,2] - shear_data_2[:,3])

        count_2 = numpy.histogram2d(mg2_1, mg2_2, bins)[0]
        total_count_2[i] = total_count_2[i] + count_2

for i in my_task_list:
    chisq_arr[0, i] = get_chisq(total_count_1[i], mg_bin_num)
    chisq_arr[1, i] = get_chisq(total_count_2[i], mg_bin_num)

h5f_1.close()
h5f_2.close()

comm.Barrier()

t2 = time.time()
if rank == 0:
    numpy.savez(data_path + "/chisq_%d_%d_%s_narrow.npz"%(tag_1, tag_2, data_type), chi_guess_bin, chisq_arr)
    print("%d %d %.2f"%(tag_1, tag_2, t2-t1))
comm.Barrier()





