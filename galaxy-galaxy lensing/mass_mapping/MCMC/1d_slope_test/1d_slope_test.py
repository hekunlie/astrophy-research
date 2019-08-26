import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append("E:/Github/astrophy-research/mylib")
path.append('%s/work/mylib/'%my_home)
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import numpy
import matplotlib.pyplot as plt
from matplotlib import cm
from mpi4py import MPI
import tool_box


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


def plot_3d(fig, row, col, seq, x, y, z, labels):
    label_size = 20
    ax = fig.add_subplot(row, col, seq, projection='3d')
    surf = ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=cm.coolwarm,linewidth=0, antialiased=False)
    ax.set_xlabel(labels[0], fontsize=label_size)
    ax.set_ylabel(labels[1], fontsize=label_size)
    ax.set_zlabel(labels[2], fontsize=label_size)
    fig.colorbar(surf, shrink=0.5, aspect=10)

def find_patch(arr,thresh):
    ny, nx = arr.shape
    x1, x2 = -1, -1
    y1, y2 = -1, -1
    for i in range(int(ny/2)):
        if arr[i].min() <= thresh and y1 < 0:
            y1 = i
        if arr[ny-i-1].min() <= thresh and y2 < 0:
            y2 = ny-i
    for i in range(int(nx/2)):
        if arr[:,i].min() <= thresh and x1 < 0:
            x1 = i
        if arr[:, nx-i-1].min() <= thresh and x2 < 0:
            x2 = nx-i
    print(y1, y2, x1, x2)

    if min(x1, x2, y1, y2) >= 0 and x2-x1 > 1 and y2-y1>1:
        return y1, y2, x1, x2
    else:
        return 0, ny, 0, nx



seed = 1214
total_num = 250000
uni_num = 150000

chisq_nm = "chisq_asy.png"
PDF_nm = "PDF_asy.png"
cache_nm = "chisq_asy.npz"

fq = Fourier_Quad(6,12)
rng = numpy.random.RandomState(seed)

bin_num = 8
bin_num2 = int(bin_num/2)
inverse = range(int(bin_num2 - 1), -1, -1)

bin_num_test = 8
bin_num2_test = int(bin_num_test/2)
inverse_test = range(int(bin_num2_test - 1), -1, -1)

hist_bin = 20

nx, ny = 71, 71
a1_range = numpy.linspace(-0.1, 0.1, nx)
a2_range = numpy.linspace(-0.1, 0.1, ny)
# x,y in the phase space
a1_hat, a2_hat = numpy.meshgrid(a1_range, a2_range)

itemsize = MPI.DOUBLE.Get_size()
element_num = nx*ny
if rank == 0:
    nbytes = element_num*itemsize
    nbytes_gal = total_num*itemsize
else:
    nbytes = 0
    nbytes_gal = 0
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
win2 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
win3 = MPI.Win.Allocate_shared(nbytes_gal, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
buf2, itemsize = win2.Shared_query(0)
buf3, itemsize = win3.Shared_query(0)
chisq = numpy.ndarray(buffer=buf1, dtype='d', shape=(ny,nx)) # array filled with zero
chisq_new = numpy.ndarray(buffer=buf2, dtype='d', shape=(ny,nx)) # array filled with zero
mg = numpy.ndarray(buffer=buf3, dtype='d', shape=(total_num,)) # array filled with zero

comm.Barrier()


# x = numpy.random.uniform(-8, 8, gal_num)
x_asy = rng.uniform(0, 8, int(total_num - uni_num))
x_uni = rng.uniform(-2, 8, uni_num)
x = numpy.zeros((total_num,))
x[:uni_num] = x_uni
x[uni_num:] = x_asy

a1, a2 = 0, -0.01
shear_slope_1d = a1 + a2*x
print(shear_slope_1d.min(), shear_slope_1d.max())

if rank == 0:
    # simulation
    mg_ori = numpy.zeros((total_num,))
    for i in range(total_num):
        mg_ori[i] = rng.normal(0, 0.3, 1)
        mg[i] = rng.normal(shear_slope_1d[i], 0.3, 1)
        # mg[i] = mg_ori[i] + shear_slope_1d[i]

comm.Barrier()

mnu = numpy.ones_like(mg)
mg_bins = fq.set_bin(mg, bin_num, 10)
mg_bins_test = fq.set_bin(mg, bin_num_test, 10)

num_ini = numpy.histogram(mg, mg_bins)[0]
n1 = num_ini[0:bin_num2][inverse]
n2 = num_ini[bin_num2:]
num_exp = (n1 + n2) / 2

task_list = [i for i in range(nx*ny)]
sub_task = tool_box.allot(task_list, cpus)[rank]

for i in sub_task:
    m, n = divmod(i, nx)
    gh = a1_range[n] + a2_range[m]*x
    # [a2,a1]-axis
    chisq[m,n] = fq.get_chisq(mg, mnu, gh, mg_bins, bin_num2,inverse, 0)
    chisq_new[m,n] = fq.get_chisq_new(mg, mnu, gh, mg_bins, bin_num2, inverse, 0, num_exp)

comm.Barrier()

if rank == 0:
    numpy.savez(cache_nm, chisq, chisq_new)

    img = Image_Plot(fig_x=12, fig_y=9)
    img.subplots(1, 1)
    chi_temp = fq.get_chisq(mg, mnu, 0, mg_bins, bin_num2, inverse, 0)
    img.axs[0][0].hist(mg, hist_bin, histtype='step',linewidth=img.plt_line_width+2, label="Before PDF, $\chi^2$=%.2f"%chi_temp)

    # num_null = numpy.histogram(mg, mg_bins_test)[0]
    # n1 = num_null[0:bin_num2_test][inverse_test]
    # n2 = num_null[bin_num2_test:]
    # chi_null = (n1-n2)/numpy.abs(n1-n2) * numpy.sum((n1-n2)**2/(n1+n2))*0.5
    # img.axs[0][1].scatter(range(bin_num2_test), chi_null,marker="p", s=140, label="Before PDF, $\chi^2$=%.2f"%chi_temp)

    a2_test = numpy.linspace(a2-0.1, a2+0.1, 5)

    for i in range(5):
        gh = a2_test[i]*x
        chi_temp = fq.get_chisq(mg, mnu, gh, mg_bins, bin_num2, inverse, 0)
        img.axs[0][0].hist(mg-gh, hist_bin, histtype='step',linewidth=img.plt_line_width+2, label="$a_2$=%.2f, $\chi^2$=%.2f"%(a2_test[i],chi_temp))

        # num_test = numpy.histogram(mg-gh, mg_bins_test)[0]
        # n1 = num_test[0:bin_num2_test][inverse_test]
        # n2 = num_test[bin_num2_test:]
        # chi_test = (n1-n2)/numpy.abs(n1-n2) * numpy.sum((n1-n2)**2/(n1+n2))*0.5
        # img.axs[0][1].scatter(range(bin_num2_test), chi_test,s=120, label="$a_2$=%.2f, $\chi^2$=%.2f"%(a2_test[i],chi_temp))
        # img.axs[0][2].plot(range(bin_num2_test), chi_test,linewidth=img.plt_line_width, label="$\hat g$=%.2f, $\chi^2$=%f"%(a2_test[i],chi_temp))
    ys = img.axs[0][0].set_ylim()
    img.axs[0][0].plot([0,0],[ys[0], ys[1]], linestyle="--", c="grey")
    img.axs[0][0].set_ylim(ys)
    img.axs[0][0].legend(fontsize=img.legend_size-2)
    img.set_label(0,0,0,"N")
    img.set_label(0,0,1,"G")
    img.axs[0][0].set_title("Fix $a_1$=0 (True $a_2$ = %.2f)"%a2,fontsize=img.xy_lb_size)
    img.save_img(PDF_nm)
    img.show_img()


    y1,y2,x1, x2 = find_patch(chisq, 20)
    y1_n,y2_n,x1_n, x2_n = find_patch(chisq_new, 20)
    print(chisq.min(), chisq_new.min())
    print(chisq.max(), chisq_new.max())

    chisq_project = numpy.zeros_like(chisq)
    chisq_new_project = numpy.zeros_like(chisq)

    marker_size = 10
    show_chisq_value = 15

    img = Image_Plot()
    img.subplots(2,2)

    norm = plt.Normalize(vmin=numpy.min(chisq), vmax=numpy.max(chisq))
    norm_new = plt.Normalize(vmin=numpy.min(chisq_new), vmax=numpy.max(chisq_new))
    norm_show = plt.Normalize(vmin=0, vmax=show_chisq_value)
    cmap = plt.get_cmap('jet')

    for i in range(ny):
        for j in range(nx):
            # plot \chi (original) squared
            cl = cmap(norm(chisq[i, j]))
            img.axs[0][0].scatter(a1_range[i], a2_range[j], marker="s",color=cl, s=marker_size)
            # zoom in
            if chisq[i,j] <= show_chisq_value:
                chisq_project[i, j] = chisq[i, j]
                cl = cmap(norm_show(chisq[i, j]))
                img.axs[1][0].scatter(a1_range[i], a2_range[j], marker="s", color=cl, s=marker_size)

            # plot \chi (original) squared
            cl = cmap(norm_new(chisq_new[i,j]))
            img.axs[0][1].scatter(a1_range[i], a2_range[j],  marker="s",color=cl, s=marker_size)
            # zoom in
            if chisq_new[i,j] <= show_chisq_value:
                chisq_new_project[i, j] = chisq_new[i, j]
                cl = cmap(norm_show(chisq_new[i, j]))
                img.axs[1][1].scatter(a1_range[i], a2_range[j], marker="s", color=cl, s=marker_size)

    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm_new = plt.cm.ScalarMappable(cmap=cmap, norm=norm_new)
    sm_show = plt.cm.ScalarMappable(cmap=cmap, norm=norm_show)
    sm._A = []
    sm_new._A = []
    sm_show._A = []
    plt.colorbar(sm,ax=img.axs[0][0])
    plt.colorbar(sm_new,ax=img.axs[0][1])

    plt.colorbar(sm_show,ax=img.axs[1][0])
    plt.colorbar(sm_show,ax=img.axs[1][1])

    for j in range(2):
        xs = img.axs[0][j].set_xlim()
        ys = img.axs[0][j].set_ylim()
        img.axs[1][j].set_xlim(xs)
        img.axs[1][j].set_ylim(ys)

        img.axs[1][j].plot([a1,a1],[ys[0], ys[1]], c="C4", linewidth=1)
        img.axs[1][j].plot([xs[0], xs[1]],[a2, a2], c="C4", linewidth=1)

        for i in range(2):
            img.set_label(0,j,i,"$a_%d$"%(2-i))
            img.set_label(1,j,i,"$a_%d$"%(2-i))
    img.save_img(chisq_nm)
    img.show_img()

    img = Image_Plot()
    img.subplots(2,2)




# fig = plt.figure(figsize=(16, 8))
# plot_3d(fig, 1,2,1,a1_hat, a2_hat, chisq,["$a_1$","$a_2$","$\chi^2$"])
# plot_3d(fig, 1,2,2,a1_hat, a2_hat, chisq_new,["$a_1$","$a_2$","$\chi^2_{new}$"])
# plt.show()

# fig = plt.figure(figsize=(16, 8))
# plot_3d(fig, 1,2,1,a1_hat[y1:y2,x1:x2], a2_hat[y1:y2,x1:x2], chisq[y1:y2,x1:x2],["$a_1$","$a_2$","$\chi^2$"])
# plot_3d(fig, 1,2,2,a1_hat[y1_n:y2_n,x1_n:x2_n], a2_hat[y1_n:y2_n,x1_n:x2_n], chisq_new[y1_n:y2_n,x1_n:x2_n],
#         ["$a_1$","$a_2$","$\chi^2_{new}$"])
# plt.show()