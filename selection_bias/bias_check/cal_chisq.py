import h5py
from mpi4py import MPI
from sys import path,argv
path.append("/home/hklee/work/mylib")
import tool_box
from Fourier_Quad import Fourier_Quad
import numpy
from plot_tool import Image_Plot


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

fq = Fourier_Quad(12,124)

shear_num = 20
n, m = divmod(shear_num, numprocs)
tasks = [i for i in range(shear_num)]

my_task = tool_box.allot(tasks,numprocs)[rank]

print(rank, my_task)

total_path = argv[1]
# # columns: [col_st, col_ed]
# col_st, col_ed = int(argv[2]), int(argv[3])
# # name of the sub-data file
# sep_data_nm = argv[4]

dst_nms = ["gauss_noise_1", "gauss_noise_2", "gauss_noise_residual",
           "moffat_noise_1", "moffat_noise_2", "moffat_noise_residual"]

guess_num = 31
g_range = [[-0.0005,0.0005],[-0.0005,0.0005],[-0.1,0.1],
           [-0.0005,0.0005],[-0.0005,0.0005],[-0.1,0.1]]

sub_num = 10
result = numpy.zeros((4, guess_num))
sub_result = numpy.zeros((sub_num*2+1, guess_num))

for ig in my_task:
    for i in range(len(dst_nms)):
        src_path = total_path + "/data_%d_%s.hdf5"%(ig,dst_nms[i])
        h5f = h5py.File(src_path,"r")
        data = h5f["/data"].value
        h5f.close()

        mg1 = data[:, 0]
        mg2 = data[:, 1]
        mn = data[:, 2]
        mu = data[:, 3]
        mnu1 = mn + mu
        mnu2 = mn - mu

        shear_guess = numpy.linspace(g_range[i][0],g_range[i][1], guess_num)
        g1_hat, g1_chisq = fq.get_chisq_range(mg1, mnu1, 8, shear_guess)
        g2_hat, g2_chisq = fq.get_chisq_range(mg2, mnu2, 8, shear_guess)

        result[0] = g1_hat
        result[1] = g1_chisq
        result[2] = g2_hat
        result[3] = g2_chisq

        img = Image_Plot(fig_x=6,fig_y=4,xpad=0.25,ypad=0.25)
        img.subplots(1,1)
        img.axs[0][0].plot(g1_hat, g1_chisq,marker="o",label="g1")
        img.axs[0][0].plot(g2_hat, g2_chisq,marker="s",label="g2")
        img.axs[0][0].legend(fontsize=img.legend_size)
        img.set_label(0,0,0,"$\chi^2$",fontsize=img.xy_lb_size)
        img.set_label(0,0,1,"$g$",fontsize=img.xy_lb_size)

        img.save_img(total_path + "/result/pic/data_%d_%s_chisq.png"%(ig,dst_nms[i]))
        img.close_img()

        if i == 2 or i == 5:

            img = Image_Plot(fig_x=6, fig_y=4, xpad=0.25, ypad=0.25)
            img.subplots(1, 2)

            total_row = mg1.shape[0]
            sub_row = divmod(total_row, sub_num)[0]
            sub_result[2*sub_num] = shear_guess

            for j in range(sub_num):
                g1_hat, g1_chisq = fq.get_chisq_range(mg1[j*sub_row:(j+1)*sub_row], mnu1[j*sub_row:(j+1)*sub_row], 8, shear_guess)
                g2_hat, g2_chisq = fq.get_chisq_range(mg2[j*sub_row:(j+1)*sub_row], mnu2[j*sub_row:(j+1)*sub_row], 8, shear_guess)

                sub_result[j] = g1_chisq
                sub_result[j+sub_num] = g2_chisq

                img.axs[0][0].plot(g1_hat, g1_chisq, marker="o")
                img.axs[0][1].plot(g2_hat, g2_chisq, marker="o")

            for j in range(2):
                img.set_label(0,j,0,"$\chi^2$",fontsize=img.xy_lb_size)
                img.set_label(0,j,1,"$g$",fontsize=img.xy_lb_size)

            img.save_img(total_path + "/result/pic/data_%d_%s_chisq_sub_sample.png" % (ig, dst_nms[i]))
            img.close_img()

        dst_path = total_path + "/result/data_%d_%s_chisq.hdf5"%(ig,dst_nms[i])
        h5f = h5py.File(dst_path,"w")
        h5f["/data"] = result
        h5f["/sub_data"] = sub_result
        h5f.close()


