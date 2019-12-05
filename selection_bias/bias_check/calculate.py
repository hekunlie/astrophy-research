import numpy
import h5py
from sys import argv, path
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
import tool_box
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

cmd = argv[1]
shear_num = int(argv[2])
total_path = argv[3]
noise_type = argv[4]

data_type = "/data_%d_"+"%s.hdf5"%noise_type
result_name = total_path + "/result/result_%s.hdf5"%noise_type
img_check_name = total_path + "/result/data_%d_"+"%s.png"%noise_type
result_img_name = total_path + "/result/result_%s.png"%noise_type


def mc_plot(mean_result,mean_mc,sym_result,sym_mc,img_name):

    img = Image_Plot(fig_x=6, fig_y=4, xpad=0.25, ypad=0.2)
    img.subplots(2, 2)

    scale = 1000

    img.axs[0][0].errorbar(mean_result[0], mean_result[1], mean_result[2],
                           marker="s", fmt=" ", label="g1", capsize=img.cap_size)
    img.axs[0][0].errorbar(mean_result[3], mean_result[4], mean_result[5],
                           marker="s", fmt=" ", label="g2", capsize=img.cap_size)
    img.axs[1][0].errorbar(mean_result[0], scale*(mean_result[1]-mean_result[0]), scale*mean_result[2],
                           marker="s", fmt=" ", label="g1", capsize=img.cap_size)
    img.axs[1][0].errorbar(mean_result[3], scale*(mean_result[4]-mean_result[3]), scale*mean_result[5],
                           marker="s", fmt=" ", label="g2", capsize=img.cap_size)

    mc_str = "Mean\n$m_1 = %.5f (%.5f)$\n$c_1 = %.5f (%.5f)$\n$m_2 = %.5f (%.5f)$\n$c_2 = %.5f (%.5f)$" \
             % (mean_mc[0, 0], mean_mc[0, 1], mean_mc[0, 2], mean_mc[0, 3],
                mean_mc[1, 0], mean_mc[1, 1], mean_mc[1, 2], mean_mc[1, 3])
    img.axs_text(0, 0, 0.8, 0.05, mc_str, text_fontsize=img.legend_size)


    img.axs[0][1].errorbar(sym_result[0], sym_result[1], sym_result[2],
                           marker="s", fmt=" ", label="g1", capsize=img.cap_size)
    img.axs[0][1].errorbar(sym_result[3], sym_result[4], sym_result[5],
                           marker="s", fmt=" ", label="g2", capsize=img.cap_size)

    img.axs[1][1].errorbar(sym_result[0], scale*(sym_result[1]-sym_result[0]), scale*sym_result[2],
                           marker="s", fmt=" ", label="g1", capsize=img.cap_size)
    img.axs[1][1].errorbar(sym_result[3], scale*(sym_result[4]-sym_result[3]), scale*sym_result[5],
                           marker="s", fmt=" ", label="g2", capsize=img.cap_size)

    mc_str = "PSF_SYM\n$m_1 = %.5f (%.5f)$\n$c_1 = %.5f (%.5f)$\n$m_2 = %.5f (%.5f)$\n$c_2 = %.5f (%.5f)$" \
             % (sym_mc[0, 0], sym_mc[0, 1], sym_mc[0, 2], sym_mc[0, 3],
                sym_mc[1, 0], sym_mc[1, 1], sym_mc[1, 2], sym_mc[1, 3])
    img.axs_text(0, 1, 0.8, 0.05, mc_str, text_fontsize=img.legend_size)


    for i in range(2):
        img.set_label(0, i, 0, "EST g")
        img.set_label(1, i, 0, "$10^3 (\hat{g} - g_{true})$")
        for j in range(2):
            img.set_label(i, j, 1, "TRUE g")
            img.axs[i][j].legend(fontsize=img.legend_size, loc="lower right")
            img.axs[i][j].plot([-0.1, 0.1], [-0.1, 0.1], ls="dashed", c="grey", alpha=0.5)
            img.axs[0][j].set_ylim(-0.05, 0.05)
            img.axs[1][j].set_ylim(-0.5, 0.5)
            img.axs[i][j].set_xlim(-0.05, 0.05)

    img.save_img(img_name)


if cmd == "calculate":
    fq = Fourier_Quad(12,124)

    itemsize = MPI.DOUBLE.Get_size()

    if rank == 0:
        nbytes1 = shear_num*6*itemsize
        nbytes2 = 2*4*itemsize
    else:
        nbytes1 = 0
        nbytes2 = 0

    win1 = MPI.Win.Allocate_shared(nbytes1, itemsize, comm=comm)
    win2 = MPI.Win.Allocate_shared(nbytes2, itemsize, comm=comm)
    win3 = MPI.Win.Allocate_shared(nbytes1, itemsize, comm=comm)
    win4 = MPI.Win.Allocate_shared(nbytes2, itemsize, comm=comm)

    buf1, itemsize = win1.Shared_query(0)
    buf2, itemsize = win2.Shared_query(0)
    buf3, itemsize = win3.Shared_query(0)
    buf4, itemsize = win4.Shared_query(0)

    # true g1, g1, g1_sig, true g2, g2, g2_sig
    mean_result_array = numpy.ndarray(buffer=buf1, dtype='d', shape=(6,shear_num)) # array filled with zero
    # [m1, m1_sig, c1, c1_sig],[m2, m2_sig, c2, c2_sig]
    mean_mc_array = numpy.ndarray(buffer=buf2, dtype='d', shape=(2, 4))

    sym_result_array = numpy.ndarray(buffer=buf3, dtype='d', shape=(6,shear_num)) # array filled with zer
    sym_mc_array = numpy.ndarray(buffer=buf4, dtype='d', shape=(2, 4))

    tasks = [i for i in range(shear_num)]
    my_task = tool_box.allot(tasks,numprocs)[rank]

    h5f = h5py.File(total_path + "/shear.hdf5", "r")
    g1 = h5f["/g1"].value
    g2 = h5f["/g2"].value
    h5f.close()

    if rank == 0:

        mean_result_array[0] = g1
        mean_result_array[3] = g2
        sym_result_array[0] = g1
        sym_result_array[3] = g2

    comm.Barrier()

    for ig in my_task:
        src_path = total_path + data_type%ig
        h5f = h5py.File(src_path, "r")
        data_ig = h5f["/data"].value
        h5f.close()

        mg1 = data_ig[:, 0]
        mg2 = data_ig[:, 1]
        mn = data_ig[:, 2]
        mu = data_ig[:, 3]
        mnu1 = data_ig[:, 2] + data_ig[:, 3]
        mnu2 = data_ig[:, 2] - data_ig[:, 3]
        print("Doing %d, %s"%(ig, src_path), data_ig.shape)

        img_check = Image_Plot()
        img_check.subplots(1,2)

        left, right = -0.06, 0.06
        dg = 0.02
        g1_mean, g1_sig_mean = fq.find_shear_mean(mg1, mn)
        g2_mean, g2_sig_mean = fq.find_shear_mean(mg2, mn)

        g1_sym, g1_sig_sym = fq.find_shear(mg1, mnu1, 8, left=left, right=right, fit_num=20, chi_gap=30,
                                           fig_ax=img_check.axs[0][0])[:2]
        g2_sym, g2_sig_sym = fq.find_shear(mg2, mnu2, 8, left=left, right=right, fit_num=20, chi_gap=30,
                                           fig_ax=img_check.axs[0][1])[:2]

        img_check.save_img(img_check_name%ig)
        img_check.close_img()

        mean_result_array[1,ig] = g1_mean
        mean_result_array[2,ig] = g1_sig_mean
        mean_result_array[4,ig] = g2_mean
        mean_result_array[5,ig] = g2_sig_mean

        sym_result_array[1,ig] = g1_sym
        sym_result_array[2,ig] = g1_sig_sym
        sym_result_array[4,ig] = g2_sym
        sym_result_array[5,ig] = g2_sig_sym

    comm.Barrier()

    if rank == 0:

        mc1 = tool_box.data_fit_numpy(mean_result_array[0], mean_result_array[1], mean_result_array[2])
        mean_mc_array[0] = mc1[0]-1, mc1[1], mc1[2], mc1[3]
        mc2 = tool_box.data_fit_numpy(mean_result_array[3], mean_result_array[4], mean_result_array[5])
        mean_mc_array[1] = mc2[0]-1, mc2[1], mc2[2], mc2[3]

        mc1 = tool_box.data_fit_numpy(sym_result_array[0], sym_result_array[1], sym_result_array[2])
        sym_mc_array[0] = mc1[0]-1, mc1[1], mc1[2], mc1[3]
        mc2 = tool_box.data_fit_numpy(sym_result_array[3], sym_result_array[4], sym_result_array[5])
        sym_mc_array[1] = mc2[0]-1, mc2[1], mc2[2], mc2[3]

        h5f = h5py.File(result_name,"w")
        h5f["/mean_result"] = mean_result_array
        h5f["/mean_mc"] = mean_mc_array
        h5f["/sym_result"] = sym_result_array
        h5f["/sym_mc"] = sym_mc_array
        h5f.close()

        mc_plot(mean_result_array, mean_mc_array,sym_result_array, sym_mc_array, result_img_name)

        print("Mean:")
        print(mean_mc_array)

        print("SYM:")
        print(sym_mc_array)

else:
    if rank == 0:

        h5f = h5py.File(result_name, "r")
        mean_result_array = h5f["/mean_result"].value
        mean_mc_array = h5f["/mean_mc"].value
        sym_result_array = h5f["/sym_result"].value
        sym_mc_array = h5f["/sym_mc"].value
        h5f.close()

        mc_plot(mean_result_array, mean_mc_array,sym_result_array, sym_mc_array, result_img_name)

        print("Mean:")
        print(mean_mc_array)

        print("SYM:")
        print(sym_mc_array)
    comm.Barrier()