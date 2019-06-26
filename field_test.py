import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import matplotlib.pyplot as plt
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
from Fourier_Quad import Fourier_Quad
import h5py
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

t1 = time.time()

field_name = argv[1]
cmd = argv[2]

data_path = "/mnt/ddnfs/data_users/hkli/CFHT/catalog/fourier_cata/%s/result_ext/%s_shear.dat"%(field_name,field_name)
# data_path = "/home/hkli/work/CFHT_from_Jun/%s_shear.cat"%field_name
data = numpy.loadtxt(data_path)

envs_path = "%s/work/envs/envs.dat"%my_home
gets_item = [["fresh_para_idx", "nstar", "0"], ["fresh_para_idx", "flux_alt", "0"],
             ["fresh_para_idx", "ra", "0"], ["fresh_para_idx", "dec", "0"],
             ["fresh_para_idx", "gf1", "0"], ["fresh_para_idx", "gf2", "0"],
             ["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
             ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
             ["fresh_para_idx", "h2", "0"], ["fresh_para_idx", "total_area", "0"]]
gets = ["get" for i in range(len(gets_item))]
para_items = tool_box.config(envs_path, gets, gets_item)

nstar_lb = int(para_items[0])
flux_alt_lb = int(para_items[1])
total_area_lb = int(para_items[11])

ra_lb = int(para_items[2])
dec_lb = int(para_items[3])

field_g1_lb = int(para_items[4])
field_g2_lb = int(para_items[5])

mg1_lb = int(para_items[6])
mg2_lb = int(para_items[7])
mn_lb = int(para_items[8])
mu_lb = int(para_items[9])
mv_lb = int(para_items[10])

# field distortion
fg1 = data[:,field_g1_lb]
fg2 = data[:,field_g2_lb]

mg1 = data[:, mg1_lb] / 2
mg2 = data[:, mg2_lb] / 2
mn = data[:, mn_lb] / 2
mu = - data[:, mu_lb] / 2
mv = - data[:, mv_lb] / 2
mnu1 = mn + mu
mnu2 = mn - mu

# nstar_lb = 11
# flux_alt_lb = 10
#
# ra_lb = 0
# dec_lb = 1
#
# field_g1_lb = 0
# field_g2_lb = 1
#
# mg1_lb = 3
# mg2_lb = 4
# mn_lb = 5
# mu_lb = 6
# mv_lb = 7
#
#
# mg1 = data[:, mg1_lb] / 2
# mg2 = data[:, mg2_lb] / 2
# mn = data[:, mn_lb] / 2
# mu = - data[:, mu_lb] / 2
# mv = - data[:, mv_lb] / 2
# mnu1 = mn + mu
# mnu2 = mn - mu
#
# # field distortion
# fg1 = numpy.zeros_like(mg1)
# fg2 = numpy.zeros_like(mg1)

fq = Fourier_Quad(12, 4)

# plot g1 vs field_g1 & g2 vs field_g2
if cmd == "field":
    fg1_bin_num = 40
    fg2_bin_num = 60
    fg1_bin = numpy.linspace(-0.005, 0.005, fg1_bin_num+1)
    fg2_bin = numpy.linspace(-0.0075, 0.0075, fg2_bin_num+1)

    # cutoff
    idx1 = numpy.abs(fg1) <= 0.005
    idx2 = numpy.abs(fg2) <= 0.0075
    idx3 = data[:,flux_alt_lb] >= 3.5

    idx = idx1 & idx2 & idx3

    if rank == 0:
        print(mg1_lb, mg2_lb, mn_lb, mu_lb, field_g1_lb, field_g2_lb)
        print("Totally %d galaxies. %d have been selected"%(idx.shape[0],idx.sum()))
    # shear estimators
    fg1 = fg1[idx]
    fg2 = fg2[idx]
    mg1 = mg1[idx]
    mg2 = mg2[idx]
    mn = mn[idx]
    mu = mu[idx]
    mnu1 = mnu1[idx]
    mnu2 = mnu2[idx]

    itemsize = MPI.DOUBLE.Get_size()
    if rank == 0:
        # bytes for 10 double elements
        nbytes = fg2_bin_num*4*itemsize
    else:
        nbytes = 0

    win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
    buf1, itemsize = win1.Shared_query(0)
    result = numpy.ndarray(buffer=buf1, dtype='d', shape=(4, fg2_bin_num))
    comm.Barrier()

    my_g1 = tool_box.allot([i for i in range(fg1_bin_num)], cpus)[rank]
    for i in my_g1:
        idx_1 = fg1 >= fg1_bin[i]
        idx_2 = fg1 < fg1_bin[i+1]
        idx_t = idx_1 & idx_2
        g1, g1_sig = fq.fmin_g_new(mg1[idx_t], mnu1[idx_t], 8, fit_num=20)[:2]
        result[0,i] = g1
        result[1,i] = g1_sig

    my_g2 = tool_box.allot([i for i in range(fg2_bin_num)], cpus)[rank]
    for i in range(fg2_bin_num):
        idx_1 = fg2 >= fg2_bin[i]
        idx_2 = fg2 < fg2_bin[i+1]
        idx_t = idx_1 & idx_2
        g2, g2_sig = fq.fmin_g_new(mg2[idx_t], mnu2[idx_t], 8, fit_num=20)[:2]
        result[2,i] = g2
        result[3,i] = g2_sig

    comm.Barrier()

    t2 = time.time()

    if rank == 0:
        img = Image_Plot()
        img.subplots(1,2)

        img.axs[0][0].errorbar(fg1_bin[:fg1_bin_num], result[0, :fg1_bin_num], result[1, :fg1_bin_num], capsize=3, c="C1", label="g1")
        img.axs[0][0].scatter(fg1_bin[:fg1_bin_num], result[0, :fg1_bin_num], c="C1", marker="s")

        xs = img.axs[0][0].set_xlim()
        ys = img.axs[0][0].set_ylim()
        max_x = max([abs(xs[0]), abs(xs[1])])
        max_y = max([abs(ys[0]), abs(ys[1])])
        max_range = max(max_x, max_y)

        img.axs[0][0].plot([-max_range, max_range], [-max_range, max_range], linestyle="--", c="grey")
        img.axs[0][0].legend(loc="upper left")
        img.axs[0][0].set_xlim(-max_x, max_x)
        img.axs[0][0].set_ylim(-max_y, max_y)

        img.set_label(0, 0, 0, "$g_1$")
        img.set_label(0, 0, 1, "True $g_1$")

        img.axs[0][1].errorbar(fg2_bin[:fg2_bin_num], result[2, :fg2_bin_num], result[3, :fg2_bin_num], capsize=3, c="C1", label="g2")
        img.axs[0][1].scatter(fg2_bin[:fg2_bin_num], result[2, :fg2_bin_num], c="C1", marker="s")

        xs = img.axs[0][1].set_xlim()
        ys = img.axs[0][1].set_ylim()
        max_x = max([abs(xs[0]), abs(xs[1])])
        max_y = max([abs(ys[0]), abs(ys[1])])
        max_range = max(max_x, max_y)

        img.axs[0][1].plot([-max_range, max_range], [-max_range, max_range], linestyle="--", c="grey")
        img.axs[0][1].legend(loc="upper left")
        img.axs[0][1].set_xlim(-max_x, max_x)
        img.axs[0][1].set_ylim(-max_y, max_y)

        img.set_label(0, 1, 0, "$g_2$")
        img.set_label(0, 1, 1, "True $g_2$")
        img.save_img("%s.png" % field_name)
        print("%.2f"%(t2-t1))

# plot shear map
if cmd == "map":
    t1 = time.time()

    # correct the field distortion
    mg1 = mg1 - fg1 * (mn + mu) - fg2 * mv
    mg2 = mg2 - fg2 * (mn - mu) - fg1 * mv

    ra = data[:, ra_lb]
    dec = data[:, dec_lb]

    ra_min, ra_max = ra.min(), ra.max()
    dec_min, dec_max = dec.min(), dec.max()

    # set up bins for RA and DEC
    ny, nx = int(argv[3]), int(argv[4])
    ra_bin = numpy.linspace(ra_min, ra_max, nx+1)
    dec_bin = numpy.linspace(dec_min, dec_max, ny+1)
    grid_num = ny*nx

    itemsize = MPI.DOUBLE.Get_size()
    if rank == 0:
        # bytes for 10 double elements
        nbytes = grid_num*5*itemsize
    else:
        nbytes = 0

    win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
    buf1, itemsize = win1.Shared_query(0)
    # [[g1], [g2]]
    shear = numpy.ndarray(buffer=buf1, dtype='d', shape=(ny, 5*nx))
    comm.Barrier()

    my_grid = tool_box.allot([(i, j) for i in range(ny) for j in range(nx)], cpus)[rank]

    for grid_id in my_grid:
        iy, ix = grid_id

        idx_1 = dec >= dec_bin[iy]
        idx_2 = dec < dec_bin[iy+1]
        idx_3 = ra >= ra_bin[ix]
        idx_4 = ra < ra_bin[ix+1]

        idx = idx_1 & idx_2 & idx_3 & idx_4

        if idx.sum() > 0:
            g1, g1_sig = fq.fmin_g_new(mg1[idx], mnu1[idx], 8, fit_num=20)[:2]
            g2, g2_sig = fq.fmin_g_new(mg2[idx], mnu2[idx], 8, fit_num=20)[:2]

            shear[iy, ix] = g1
            shear[iy, ix + nx] = g1_sig
            shear[iy, ix + 2*nx] = g2
            shear[iy, ix + 3*nx] = g2_sig
            shear[iy, ix + 4*nx] = idx.sum()

    comm.Barrier()
    t2 = time.time()
    if rank == 0:
        img_g = Image_Plot(fig_x=8, fig_y=8)
        img_g.subplots(2, 3)

        fig_g1 = img_g.axs[0][0].imshow(shear[:, :nx])
        plt.colorbar(fig_g1, ax=img_g.axs[0][0])

        fig_g2 = img_g.axs[0][1].imshow(shear[:, 1*nx:2*nx])
        plt.colorbar(fig_g2, ax=img_g.axs[0][1])

        print(shear[:, 1*nx:2*nx]*numpy.sqrt(shear[:, 4*nx:5*nx]))
        fig_g3 = img_g.axs[1][0].imshow(shear[:, 2*nx:3*nx])
        plt.colorbar(fig_g3, ax=img_g.axs[1][0])

        fig_g4 = img_g.axs[1][1].imshow(shear[:, 3*nx:4*nx])
        plt.colorbar(fig_g4, ax=img_g.axs[1][1])
        print(shear[:, 3*nx:4*nx] * numpy.sqrt(shear[:, 4 * nx:5 * nx]))
        fig_g4 = img_g.axs[1][2].imshow(shear[:, 4*nx:5*nx])
        plt.colorbar(fig_g4, ax=img_g.axs[1][2])


        img_g.save_img("shear_map_%s.png"%field_name)
        #
        # img_map = Image_Plot(fig_x=15, fig_y=15)
        # img_map.subplots(1, 1)
        # for i in range(ny):
        #     for j in range(nx):
        #         img
        print("%.2f"%(t2-t1))





