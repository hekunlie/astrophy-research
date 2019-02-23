import matplotlib
matplotlib.use("Agg")
import numpy
import matplotlib.pyplot as plt
import os
from sys import path
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
path.append('%s/work/fourier_quad/'%my_home)
import time
from Fourier_Quad import Fourier_Quad
from sys import argv
from mpi4py import MPI
import tool_box
import h5py


select_cri, nstar_thres, area_thres, result_source = argv[1], int(argv[2]), int(argv[3]), argv[4]

ts = time.time()

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

envs_path = "%s/work/envs/envs.dat" % my_home
fresh = "fresh_para_idx"
get_contents = [['cfht', "cfht_path_result", '1'],['cfht', "cfht_path_data", '1'],['cfht', "cfht_path_cut", '1'],
                ["%s" % fresh, "nstar", "1"], ["%s" % fresh, "total_area", "1"],
                ["%s" % fresh, "flux2", '1'], ["%s" % fresh, "flux_alt", '1']]
gets = ["get" for i in range(len(get_contents))]
path_items = tool_box.config(envs_path, gets, get_contents)

result_path, data_path, cut_path = path_items[0:3]
nstar_lb, total_area_lb = int(path_items[3]), int(path_items[4])
flux2_lb, flux_alt_lb = int(path_items[5]), int(path_items[6])

para_path = "./para.ini"
get_contents = [['field', "g1_num", '1'], ['field', "g2_num", '1'], ['field', "g1_s", '1'],
                ['field', "g1_e", '1'], ['field', "g2_s", '1'], ['field', "g2_e", '1']]
gets = ["get" for i in range(len(get_contents))]
path_items = tool_box.config(para_path, gets, get_contents)

# parameters for segment
g1_num = int(path_items[0])
g2_num = int(path_items[1])
g1_s = float(path_items[2])
g1_e = float(path_items[3])
g2_s = float(path_items[4])
g2_e = float(path_items[5])


fg1 = numpy.linspace(g1_s, g1_e, g1_num)
fg2 = numpy.linspace(g2_s, g2_e, g2_num)
dfg1 = (fg1[1] - fg1[0])/2
dfg2 = (fg2[1] - fg2[0])/2
fgs = [fg1, fg2]
dfgs = [dfg1, dfg2]
gnums = [g1_num, g2_num]
cut_num = 20

# nstar total_area flux2 flux_alt gf1 gf2 g1(G1) g2(G2) de(N) h1(U) h2(V)
# 0     1          2     3         4   5   6      7      8    9     10
select_idx = {"flux2": 2, "flux_alt": 3}

# the column in original catalog
ori_data_labels = {"flux2": flux2_lb, "flux_alt": flux_alt_lb}

fq = Fourier_Quad(48, 123)

itemsize = MPI.DOUBLE.Get_size()

if rank == 0:
    nb1 = int(itemsize*g1_num*3*cut_num)
    nb2 = int(itemsize*g2_num*3*cut_num)
    nb3 = int(itemsize*cut_num)
else:
    nb1 = 0
    nb2 = 0
    nb3 = 0
win1 = MPI.Win.Allocate_shared(nb1, itemsize, comm=comm)
win2 = MPI.Win.Allocate_shared(nb2, itemsize, comm=comm)
win3 = MPI.Win.Allocate_shared(nb3, itemsize, comm=comm)

buf1, itemsize = win1.Shared_query(0)
buf2, itemsize = win2.Shared_query(0)
buf3, itemsize = win3.Shared_query(0)

# [[measured g] -- cut_num rows
#  [sigma]      -- cut_num rows
#  [number]     -- cut_num rows   ]
# the layout of the result array in results
#       g1_1  g1_2  ... g1_i ..
# cut_1
# cut_2
# .
# .
# .
results = [numpy.ndarray(buffer=buf1, dtype='d', shape=(int(3*cut_num), g1_num)),
           numpy.ndarray(buffer=buf2, dtype='d', shape=(int(3*cut_num), g2_num))]

# rank 0 calculates the cutoff thresholds
# and saves they in temp file for all threads
# the select_scale can be seen by every threads
select_scale = numpy.ndarray(buffer=buf3, dtype='d', shape=(cut_num,))
if rank == 0:
    # read the original catalog
    ori_data_path = data_path + "cata_%s.hdf5"%result_source
    h5f = h5py.File(ori_data_path, "r")
    ori_data = h5f["/total"].value
    h5f.close()

    # choose the needed data by nstar- and total_area-cutoff
    nstar = ori_data[:,nstar_lb]
    total_area = ori_data[:, total_area_lb]
    idx_1 = nstar >= nstar_thres
    idx_2 = total_area >= area_thres

    select_data = ori_data[:, ori_data_labels[select_cri]][idx_1&idx_2]

    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111)
    ax.hist(select_data[select_data<80], 100)
    ys = ax.set_ylim()
    xs = ax.set_xlim()
    # calculate the cutoff threshold
    select_data_sort = numpy.sort(select_data)
    data_num = len(select_data)
    select_data_step = int(data_num/ cut_num)
    for i in range(cut_num):
        select_scale[i] = select_data_sort[i * select_data_step]
        # plot the cutoff threshold on the histogram
        ax.plot([select_data_sort[i * select_data_step],select_data_sort[i * select_data_step]],
                [ys[0], ys[1]], linestyle="--",c="grey",linewidth=1)
    numpy.savez(cut_path + "%s/cut_threshold.npz"%select_cri, select_scale)
    ax.set_xlim(xs[0], xs[1]+5)
    plt.savefig(cut_path + "%s/%s.pdf"%(select_cri,select_cri))
    plt.close()

    print("%7d, %8.4f, %8.4f, %7d, %8.4f, %8.4f"%(g1_num, g1_s, g1_e, g2_num, g2_s, g2_e))
    print(select_cri, select_scale)
comm.Barrier()

# start cutoff
mcs = [[], []]
# fg1 & fg2
for i in range(2):
    cat_data_path = data_path + "cata_%s_g%d.hdf5"%(result_source, i+1)

    # tasks distribution
    # if cpus <= total tasks, each core will get at least one task,
    # labeled by a number (0 ~ total_task)
    # if cpus > total tasks, some cores will get one task,
    # and the other will get nothing, a empty task list '[]'.
    tasks = tool_box.allot([i for i in range(gnums[i])], cpus)[rank]
    if len(tasks) != 0:
        for ig in tasks:
            # read the data
            set_name = "/fg%d_%d" % (i + 1, ig)
            h5f = h5py.File(cat_data_path, "r")
            data = h5f[set_name].value
            h5f.close()

            nstar = data[:, 0]
            idx_s = nstar >= nstar_thres
            area = data[:, 1]
            idx_a = area >= area_thres

            selected_data = data[idx_s & idx_a]

            if i == 0:
                # g1
                mg = selected_data[:, 6]
                nu = selected_data[:, 8] + selected_data[:, 9]
            else:
                # g2
                mg = selected_data[:, 7]
                nu = selected_data[:, 8] - selected_data[:, 9]

            for ic in range(cut_num):
                # the cutoff
                idx = selected_data[:, select_idx[select_cri]] >= select_scale[ic]

                mgs = mg[idx]
                nus = nu[idx]
                source_num = len(mgs)
                g_h, g_sig = fq.fmin_g_new(mgs, nus, bin_num=8)

                results[i][ic*3, ig] = g_h
                results[i][ic*3+1, ig] = g_sig
                results[i][ic*3+2, ig] = source_num

    comm.Barrier()
    if rank == 0:
        # measured g
        ch1 = [i * 3 for i in range(cut_num)]
        all_measured_fg = results[i][ch1]
        # sigma
        ch2 = [i * 3 + 1 for i in range(cut_num)]
        all_sigma_fg = results[i][ch2]
        # number
        ch3 = [i * 3 + 2 for i in range(cut_num)]
        all_source_num = results[i][ch3]
        x = numpy.linspace(-0.02, 0.02, 100)
        for ic in range(cut_num):
            # mc = (1+m, sig_, c, sig_c)
            mc = tool_box.data_fit(fgs[i], all_measured_fg[ic], all_sigma_fg[ic])
            mcs[i].append(mc)

            # plot the result
            fig = plt.figure(figsize=(8,8))
            ax = fig.add_subplot(111)

            lb = "fg_%d: %d \n$10^2$m: %.2f (%.2f), \n$10^4$c: %.2f (%.2f)"\
                 %(i+1,all_source_num[ic].sum(), 100*(mc[0]-1), 100*mc[1], 10**4*mc[2], 10**4*mc[3])
            ax.errorbar(fgs[i], all_measured_fg[ic], all_sigma_fg[ic], c='red',alpha=0.8, capsize=3, label=lb)
            ax.plot(mc[0]*x+mc[2], c='blue', linestyle="-.", linewidth=1)
            x1, x2 = fgs[i][0], fgs[i][-1]
            ax.plot([x1,x2],[x1,x2], c='black', linewidth=1)
            ax.set_xlim(x1 - 0.3 * (x2 - x1), x2 + 0.3 * (x2 - x1))
            ax.set_ylim(x1 - 0.3 * (x2 - x1), x2 + 0.3 * (x2 - x1))

            ax.legend()
            pic_path = cut_path + "%s/g%d_%.2f.png"%(select_cri, i+1, select_scale[ic])
            plt.savefig(pic_path)
            plt.close()

    comm.Barrier()


if rank == 0:
    mc1 = numpy.array(mcs[0]).T
    mc2 = numpy.array(mcs[1]).T
    mc_path = cut_path + "%s/total.npz"%select_cri
    numpy.savez(mc_path, mc1, mc2, results[0], results[1])
    # mc1 = numpy.load(mc_path)['arr_0']
    # mc2 = numpy.load(mc_path)['arr_1']

    x1 = 0
    x2 = 1
    x_coord = [ic * 1. / cut_num for ic in range(cut_num)]

    fig = plt.figure(figsize=(10, 5))

    ax1 = fig.add_subplot(121)
    ax1.errorbar(x_coord, mc1[0] - 1, mc1[1], c='coral', capsize=3, label='m1')
    ax1.errorbar(x_coord, mc2[0] - 1, mc2[1], c='royalblue', capsize=3, label='m2')
    ax1.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='grey')
    ax1.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax1.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax1.set_xlabel("Cutoff")
    ax1.legend()
    # ax1.set_xscale('log')
    ax1.set_ylabel("m")

    ax2 = fig.add_subplot(122)
    ax2.errorbar(x_coord, mc1[2], mc1[3], c='coral', capsize=3, label='c1')
    ax2.errorbar(x_coord, mc2[2], mc2[3], c='royalblue', capsize=3, label='c2')
    ax2.plot([x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1)], [0, 0], c='grey')
    ax2.set_xlim(x1 - 0.05 * (x2 - x1), x2 + 0.05 * (x2 - x1))
    ax2.yaxis.get_major_formatter().set_powerlimits((1, 2))
    ax2.set_xlabel("Cutoff")
    ax2.legend()
    # ax2.set_xscale('log')
    ax2.set_ylabel("c")

    namep = cut_path + "%s/total.pdf"%select_cri
    plt.savefig(namep)
    plt.close()

comm.Barrier()

te = time.time()
if rank == 0:
    print(te - ts)
