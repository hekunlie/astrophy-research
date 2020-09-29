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
# cut_name = argv[2]
sub_field_tag = int(argv[2])


fq = Fourier_Quad(12,124)
cata_name = total_path + "/data/cutoff_%d.hdf5"%sub_field_tag

h5f = h5py.File(cata_name, "r")
gf1_bin = h5f["/gf1_bin"][()]
gf2_bin = h5f["/gf2_bin"][()]
# cut_bin = h5f["/%s"%cut_name][()]
h5f.close()

shear_num1 = gf1_bin.shape[0]
shear_num2 = gf2_bin.shape[0]

# cut_num = len(cut_bin) - 1
cut_num = 1
tasks = [i for i in range(max(shear_num1, shear_num2))]

my_task = tool_box.alloc(tasks, cpus)[rank]


itemsize = MPI.FLOAT.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes1 = shear_num1 *cut_num*itemsize*2
    nbytes2 = shear_num2 *cut_num*itemsize*2
else:
    nbytes1 = 0
    nbytes2 = 0

win1 = MPI.Win.Allocate_shared(nbytes1, itemsize, comm=comm)
win2 = MPI.Win.Allocate_shared(nbytes2, itemsize, comm=comm)

buf1, itemsize = win1.Shared_query(0)
buf2, itemsize = win2.Shared_query(0)

result1 = numpy.ndarray(buffer=buf1, dtype='float32', shape=(int(2*shear_num1), cut_num))
result2 = numpy.ndarray(buffer=buf2, dtype='float32', shape=(int(2*shear_num2), cut_num))

for shear_id in my_task:
    h5f = h5py.File(cata_name, "r")

    if shear_id < shear_num1:
        mg1_all = h5f["/%d/mg_1"%shear_id][()]
        mn1_all = h5f["/%d/mn_1"%shear_id][()]
        mu1_all = h5f["/%d/mu_1"%shear_id][()]

        gh1, gh1_sig = fq.find_shear(mg1_all, mn1_all + mu1_all, 8)[:2]

        result1[shear_id, 0] = gh1
        result1[shear_num1 + shear_id, 0] = gh1_sig


        # cut_1 = h5f["/%d/%s_1"%(shear_id,cut_name)][()]
        # weight_1 = h5f["/%d/flux2_alt_1"%shear_id][()]
        #
        # for ic in range(cut_num):
        #     idx = cut_1 >= cut_bin[ic]
        #     weight = 1./weight_1[idx]
        #
        #     mg1 = mg1_all[idx]
        #     mn1 = mn1_all[idx]
        #     mu1 = mu1_all[idx]
        #
        #     # gh1, gh1_sig = fq.find_shear_mean(mg1, mn1,weight)
        #     gh1, gh1_sig = fq.find_shear(mg1, mn1 + mu1, 8)[:2]
        #
        #     result1[shear_id, ic] = gh1
        #     result1[shear_num1 + shear_id, ic] = gh1_sig
        #     print("shear: %d cut: %d. %.3f ~ %.3f. gf1: %.5f. %.5f(%.5f))"
        #           %(shear_id, ic, cut_bin[ic],cut_bin[ic+1], gf1_bin[shear_id], gh1, gh1_sig))
    if shear_id < shear_num2:
        mg2_all = h5f["/%d/mg_2"%shear_id][()]
        mn2_all = h5f["/%d/mn_2"%shear_id][()]
        mu2_all = h5f["/%d/mu_2"%shear_id][()]

        gh2, gh2_sig = fq.find_shear(mg2_all, mn2_all - mu2_all, 8)[:2]

        result2[shear_id, 0] = gh2
        result2[shear_num2 + shear_id, 0] = gh2_sig

        #
        # cut_2 = h5f["/%d/%s_2"%(shear_id,cut_name)][()]
        # weight_2 = h5f["/%d/flux2_alt_2" % shear_id][()]
        #
        # for ic in range(cut_num):
        #     idx = cut_2 >= cut_bin[ic]
        #     weight = 1./weight_2[idx]
        #
        #     mg2 = mg2_all[idx]
        #     mn2 = mn2_all[idx]
        #     mu2 = mu2_all[idx]
        #
        #     # gh2, gh2_sig = fq.find_shear_mean(mg2, mn2,weight)
        #     gh2, gh2_sig = fq.find_shear(mg2, mn2-mu2, 8)[:2]
        #
        #     result2[shear_id, ic] = gh2
        #     result2[shear_num2 + shear_id, ic] = gh2_sig
        #     print("shear: %d cut: %d. %.3f ~ %.3f. gf2: %.5f. %.5f(%.5f))"
        #           %(shear_id, ic, cut_bin[ic],cut_bin[ic+1], gf2_bin[shear_id], gh2, gh2_sig))

    h5f.close()
        # print("shear: %d cut: %d. gf1: %.5f. %.5f(%.5f) %.5f(%.5f)"%(shear_id, ic, gf_bin[shear_id], gh1, gh1_sig, gh2, gh2_sig))
comm.Barrier()

if rank == 0:

    mc1_arr = numpy.zeros((4, cut_num))
    mc2_arr = numpy.zeros((4, cut_num))

    for ic in range(cut_num):
        gh1 = result1[:shear_num1,ic]
        gh1_sig = result1[shear_num1:,ic]
        gh2 = result2[:shear_num2,ic]
        gh2_sig = result2[shear_num2:,ic]

        mc1 = tool_box.data_fit(gf1_bin, gh1, gh1_sig)
        mc1 = numpy.array(mc1)
        mc1[0] = mc1[0] - 1

        mc2 = tool_box.data_fit(gf2_bin, gh2, gh2_sig)
        mc2 = numpy.array(mc2)
        mc2[0] = mc2[0] - 1

        mc1_arr[:, ic] = mc1
        mc2_arr[:, ic] = mc2

        img = Image_Plot()
        img.subplots(1, 2)
        img.axs[0][0].errorbar(gf1_bin, result1[:shear_num1,ic], result1[shear_num1:,ic],c="C1",
                               capsize=3, marker="o", mfc="none", label="g1")
        img.axs[0][1].errorbar(gf2_bin, result2[:shear_num2,ic], result2[shear_num2:,ic],c="C1",
                               capsize=3, marker="v", mfc="none", label="g2")

        img.axs[0][0].plot([-0.01,0.01], [-0.01,0.01], c="k",ls="dashed")
        img.axs[0][1].plot([-0.01,0.01], [-0.01,0.01], c="k",ls="dashed")

        text_str = "m1: %.5f(%.5f)\nc1:%.5f(%.5f)"\
                   %(mc1[0],mc1[1],mc1[2],mc1[3])
        # img.axs_text(0,0,0.87,0.03,text_str)
        img.axs[0][0].set_title(text_str)

        text_str = "m2:%.5f(%.5f)\nc2:%.5f(%.5f)"\
                   %(mc2[0],mc2[1],mc2[2],mc2[3])
        img.axs[0][1].set_title(text_str)

        img.axs[0][0].set_xlim(-0.0065, 0.0065)
        img.axs[0][1].set_xlim(-0.0065, 0.0065)
        xlabel = ["-0.005", "-0.0025", "0", "0.0025", "0.005"]
        img.set_ticklabel_str(0, 0, 1, numpy.linspace(-0.005, 0.005, 5), xlabel)
        img.set_ticklabel_str(0, 1, 1, numpy.linspace(-0.005, 0.005, 5), xlabel)
        # img.axs[0][0].set_ylim(-0.01,0.01)
        img.axs[0][0].legend(loc="lower right")
        img.axs[0][1].legend(loc="lower right")
        # img.save_img("./result_pic/cut_%d_%s.png" %(ic,cut_name))
        img.save_img(total_path + "/result_pic/%d.png" %sub_field_tag)
        img.close_img()
    #
    numpy.savez(total_path + "/result_pic/result_%d.npz"%sub_field_tag, result1,mc1_arr,result2,mc2_arr)
    #
    # img = Image_Plot()
    # img.subplots(1,2)
    #
    # img.axs[0][0].errorbar(range(cut_num), mc1_arr[0], mc1_arr[1], capsize=3, marker="o", mfc="none", label="m1")
    # img.axs[0][0].errorbar(range(cut_num), mc2_arr[0], mc2_arr[1], capsize=3, marker="v", mfc="none", label="m2")
    # img.axs[0][0].legend()
    # img.axs[0][1].errorbar(range(cut_num), mc1_arr[2], mc1_arr[3], capsize=3, marker="o", mfc="none", label="c1")
    # img.axs[0][1].errorbar(range(cut_num), mc2_arr[2], mc2_arr[3], capsize=3, marker="v", mfc="none", label="c2")
    # img.axs[0][1].legend()
    # img.save_img("./result_pic/cut_%s.png"%cut_name)

    # print(result[:shear_num])
    # print(result[int(2*shear_num):int(3*shear_num)])
    # print(mc1_arr)
    # print(mc2_arr)

comm.Barrier()
