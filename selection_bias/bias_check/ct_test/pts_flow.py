from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
path.append("/home/hklee/work/mylib")
path.append("/home/hkli/work/mylib")
import numpy
from plot_tool import Image_Plot
import h5py
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


h5f = h5py.File("./sample_data.hdf5", "r")
shear_est_nf = h5f["/noise_free"][()]
shear_est_nr = h5f["/noise_residual"][()]
shear_est_n = h5f["/noisy"][()]
shear_est_ct = h5f["/cross_term"][()]
shear_est_ct_est = h5f["/pure_cross_term_est"][()]
h5f.close()

datas = [shear_est_nf, shear_est_nf + shear_est_ct,
         shear_est_nf + shear_est_ct_est, shear_est_nf + shear_est_ct + shear_est_ct_est,
         shear_est_nf + shear_est_ct - shear_est_ct_est,  shear_est_nf+shear_est_nr,
         shear_est_n, shear_est_ct,
         shear_est_ct_est, shear_est_ct + shear_est_ct_est,
         shear_est_ct - shear_est_ct_est, shear_est_nr]

titles = ["NF", "NF + CT", "NF + CT_est", "NF + CT + CT_est", "NF + CT - CT_est","NF + NR",
          "NY", "CT", "CT_est", "CT + CT_est","CT - CT_est", "NR"]
total_num = shear_est_ct.shape[0]
if rank == 0:
    print(total_num)

xylim = 0
xy_lock = -1
move_steps = 100
gh = numpy.linspace(0, 0.1, move_steps)

m,n = divmod(move_steps,numprocs)
st,ed = int(m*rank), int(m*(rank+1))

img_col = 6
img_row = 3

for k in range(st,ed):

    img = Image_Plot(fig_x=3, fig_y=3, xpad=0.4, ypad=0.4,xy_lb_size=12)
    img.subplots(img_row, img_col)

    for i in range(img_row):
        for j in range(img_col):

            if i < img_row - 1:
                tag = i * img_col + j
                G1_hat = datas[tag][:150, 0] + gh[k] * (datas[tag][:150, 2] + datas[tag][:150, 3])
                img.axs[i][j].scatter(G1_hat,datas[tag][:150,1],
                                      s=6, label="g1=%.3f"%gh[k])
                img.set_sci_ax(i,j,0)
                img.set_sci_ax(i,j,1)

                if xy_lock < 0:
                    ys_ = numpy.abs(img.axs[i][j].set_ylim()).max()
                    xs_ = numpy.abs(img.axs[i][j].set_xlim()).max()

                    if ys_>xylim:
                        xylim=ys_
                    if xs_>xylim:
                        xylim=xs_
            else:
                tag = (i-1) * img_col + j
                G1_hat = datas[tag][:, 0] + gh[k] * (datas[tag][:, 2] + datas[tag][:, 3])
                img.axs[i][j].hist(G1_hat, 50, label="g1=%.3f"%gh[k])
                img.set_sci_ax(i, j, 1)
                img.set_sci_ax(i, j, 0)

    xy_lock = 1
    for i in range(img_row):
        for j in range(img_col):

            if i < img_row - 1:
                tag = i * img_col + j
                img.axs[i][j].plot([-xylim,xylim],[0,0],ls="dotted",c="grey",alpha=0.5)
                img.axs[i][j].plot([0, 0], [-xylim, xylim], ls="dotted", c="grey", alpha=0.5)
                img.axs[i][j].set_ylim(-xylim, xylim)
                img.set_label(i,j,0,"G2",fontsize=img.xy_lb_size)
                img.axs[i][j].set_title(titles[tag])

            else:
                tag = (i - 1) * img_col + j
                img.set_label(i, j, 0, "Num",fontsize=img.xy_lb_size)
                img.axs[i][j].set_title(titles[tag])
                ys = img.axs[i][j].set_ylim()
                img.axs[i][j].plot([0, 0], [ys[0], ys[1]], ls="dotted", c="grey", alpha=0.5)

            img.axs[i][j].set_xlim(-xylim, xylim)
            img.set_label(i, j, 1, "G1",fontsize=img.xy_lb_size)
            img.axs[i][j].legend()

    img.save_img("./pic/move_right_%d.png"%k)
    # img.show_img()
    img.close_img()