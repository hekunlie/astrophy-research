import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import h5py
from mpi4py import MPI
import numpy
import tool_box
from Fourier_Quad import Fourier_Quad

from plot_tool import Image_Plot


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

total_path = argv[1]
shear_num = int(argv[2])
cut_name = argv[3]

fq = Fourier_Quad(12,124)
cata_name = total_path + "/cutoff.hdf5"

tasks = [i for i in range(shear_num)]

my_task = tool_box.alloc(tasks, cpus)[rank]

h5f = h5py.File(cata_name, "r")
gf_bin = h5f["/gf_bin"][()]
cut_bin = h5f["/%s"%cut_name][()]
h5f.close()
cut_num = len(cut_bin) - 1


#################### result buffer #################
# length of double
itemsize = MPI.FLOAT.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = shear_num *cut_num*itemsize*4
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf1, itemsize = win1.Shared_query(0)


# create a numpy array from buf
# code can run successfully without the following step
# buf = np.array(buf, dtype='float64', copy=False) # may be redundant
# "d" means double = 'float64'
result = numpy.ndarray(buffer=buf1, dtype='float32', shape=(int(4*shear_num), cut_num))

for shear_id in my_task:
    h5f = h5py.File(cata_name, "r")

    mg1_all = h5f["/%d/mg_1"%shear_id][()]
    mn1_all = h5f["/%d/mn_1"%shear_id][()]
    mu1_all = h5f["/%d/mu_1"%shear_id][()]
    cut_1 = h5f["/%d/%s_1"%(shear_id,cut_name)][()]
    weight_1 = h5f["/%d/flux2_alt_1"%shear_id][()]

    mg2_all = h5f["/%d/mg_2"%shear_id][()]
    mn2_all = h5f["/%d/mn_2"%shear_id][()]
    mu2_all = h5f["/%d/mu_2"%shear_id][()]
    cut_2 = h5f["/%d/%s_2"%(shear_id,cut_name)][()]
    weight_2 = h5f["/%d/flux2_alt_2" % shear_id][()]

    h5f.close()


    for ic in range(cut_num):
        idx = cut_1 >= cut_bin[ic]
        weight = weight_1[idx]
        weight[weight<0.0001] = 0.001
        weight = 1/weight
        mg1 = mg1_all[idx]
        mn1 = mn1_all[idx]
        mu1 = mu1_all[idx]
        # gh1, gh1_sig = fq.find_shear_mean(mg1, mn1, weight)
        gh1, gh1_sig = fq.find_shear(mg1, mn1, 8)[:2]

        idx = cut_2 >= cut_bin[ic]
        weight = weight_2[idx]
        weight[weight<0.0001] = 0.001
        weight = 1/weight
        mg2 = mg2_all[idx]
        mn2 = mn2_all[idx]
        mu2 = mu2_all[idx]
        # gh2, gh2_sig = fq.find_shear_mean(mg2, mn2,weight)
        gh2, gh2_sig = fq.find_shear_mean(mg2, mn2,8)[:2]


        result[shear_id,ic] = gh1
        result[shear_num + shear_id,ic] = gh1_sig
        result[int(2*shear_num) + shear_id,ic] = gh2
        result[int(3*shear_num) + shear_id,ic] = gh2_sig

        print("shear: %d cut: %d. gf1: %.5f. %.5f(%.5f) %.5f(%.5f)"%(shear_id, ic, gf_bin[shear_id], gh1, gh1_sig, gh2, gh2_sig))
comm.Barrier()

if rank == 0:


    mc1_arr = numpy.zeros((4, cut_num))
    mc2_arr = numpy.zeros((4, cut_num))

    for ic in range(cut_num):
        gh1 = result[:shear_num,ic]
        gh1_sig = result[int(shear_num):int(2*shear_num),ic]
        gh2 = result[int(2*shear_num):int(3*shear_num),ic]
        gh2_sig = result[int(3*shear_num):int(4*shear_num),ic]

        mc1 = tool_box.data_fit(gf_bin, gh1, gh1_sig)
        mc1 = numpy.array(mc1)
        mc1[0] = mc1[0] - 1

        mc2 = tool_box.data_fit(gf_bin, gh2, gh2_sig)
        mc2 = numpy.array(mc2)
        mc2[0] = mc2[0] - 1

        mc1_arr[:, ic] = mc1
        mc2_arr[:, ic] = mc2
    numpy.savez(total_path + "/result.npz", result,mc1_arr,mc2_arr)

    img = Image_Plot()
    img.subplots(1,2)

    img.axs[0][0].errorbar(range(cut_num), mc1_arr[0], mc1_arr[1], capsize=3, marker="o", mfc="none", label="m1")
    img.axs[0][0].errorbar(range(cut_num), mc2_arr[0], mc2_arr[1], capsize=3, marker="v", mfc="none", label="m2")
    img.axs[0][0].legend()
    img.axs[0][1].errorbar(range(cut_num), mc1_arr[2], mc1_arr[3], capsize=3, marker="o", mfc="none", label="c1")
    img.axs[0][1].errorbar(range(cut_num), mc2_arr[2], mc2_arr[3], capsize=3, marker="v", mfc="none", label="c2")
    img.axs[0][1].legend()
    img.save_img("cut_%s.png"%cut_name)

    # print(result[:shear_num])
    # print(result[int(2*shear_num):int(3*shear_num)])
    print(mc1_arr)
    print(mc2_arr)

comm.Barrier()
