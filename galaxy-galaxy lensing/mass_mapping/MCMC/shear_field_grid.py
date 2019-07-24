import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
import time
import matplotlib.pyplot as plt
import h5py
import MCMC_program
import tool_box
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


cmd = argv[1]
# mpirun -np 1 ....py grid
# mpirun -np 12 ....py filter cal
# mpirun -np 1 ....py filter plot
#

parent_path = "/mnt/perc/hklee/CFHT/multi_shear/cluster_field/"

envs_path = parent_path + "param_slope.dat"
contents = [['param', "RA", "1"], ['param',"DEC", "1"],
            ['param', "a1", "1"], ['param', "a2", "1"], ['param', "a3", "1"]]

var_items = tool_box.config(envs_path, ['get' for i in range(len(contents))], contents)

# the field size in unit of arcmin
delta_ra = float(var_items[0])
delta_dec = float(var_items[1])
parameters = [float(var_items[2]), float(var_items[3]), float(var_items[4])]

# read the data
fq = Fourier_Quad(10,123)
read_expo = 1
for i in range(read_expo):
    h5f = h5py.File(parent_path + "result/expo_%d_slope.hdf5"%i,"r")
    temp = h5f["/data"].value
    h5f.close()
    if i == 0:
        data = temp
    else:
        data = numpy.row_stack((data, temp))

mg1 = data[:, 0]
mg2 = data[:, 1]
mnu1 = data[:, 2] + data[:, 3]
mnu2 = data[:, 2] - data[:, 3]
shear_true = data[:, 5]
ra = data[:, 6]
dec = data[:, 7]


if cmd == "grid":
    # number of grid
    grid_nx = [5, 10]

    plt_data_coms = []
    plt_data_ms = []
    for tag, ix in enumerate(grid_nx):

        nx, ny = ix, ix
        inverse = range(nx - 1, -1, -1)

        # the grid boundary
        ra_bin = numpy.linspace(-delta_ra/2, delta_ra/2, nx+1)
        dec_bin = numpy.linspace(-delta_dec/2, delta_dec/2, ny+1)

        # the center of each grid
        grid_ra = numpy.zeros((ny, nx))
        grid_dec = numpy.zeros((ny, nx))

        # the measured shear of each grid
        shear_field_grid = numpy.zeros((ny, nx))
        shear_field_grid_err = numpy.zeros((ny, nx))

        # the mean input shear of each grid
        shear_field_grid_mean = numpy.zeros((ny, nx))

        # the galaxy number
        num_in_grid = numpy.zeros((ny, nx))

        for i in range(ny):
            idx1 = dec >= dec_bin[i]
            idx2 = dec < dec_bin[i + 1]
            cen_dec = (dec_bin[i] + dec_bin[i+1])/2
            for j in range(nx):
                idx3 = ra >= ra_bin[j]
                idx4 = ra < ra_bin[j + 1]
                cen_ra = (ra_bin[j] + ra_bin[j + 1]) / 2

                idx = idx1 & idx2 & idx3 & idx4

                num_in_grid[i, j] = idx.sum()

                gh, gh_sig = fq.fmin_g_new(mg1[idx], mnu1[idx], 8)[:2]

                shear_field_grid[i, j] = gh
                shear_field_grid_err[i, j] = gh_sig

                # center of each grid
                grid_ra[i, j] = cen_ra
                grid_dec[i, j] = cen_dec

                # mean of the input shear of each grid
                shear_field_grid_mean[i, j] = shear_true[idx].mean()

                print("%d. [%d, %d]. Num: %d. g: %.5f(%.5f).[%.5f]"%(ix, i, j,num_in_grid[i, j], gh, gh_sig, shear_field_grid_mean[i, j]))

        # the true shear on each grid center
        shear_field_grid_true = MCMC_program.shear_slope(parameters, grid_ra, grid_dec)

        # save
        data_path = parent_path + "result/grid_%d.hdf5"%ix
        h5f = h5py.File(data_path, "w")
        h5f["/ra"] = grid_ra
        h5f["/dec"] = grid_dec
        h5f["/num"] = num_in_grid
        h5f["/shear"] = shear_field_grid
        h5f["/shear_err"] = shear_field_grid_err
        h5f["/shear_mean"] = shear_field_grid_mean
        h5f["/shear_true"] = shear_field_grid_true
        h5f.close()

        # true, mean, measured, sigma, diff_percent, number
        plt_data_com = [shear_field_grid_true, shear_field_grid_mean, shear_field_grid,
                         shear_field_grid_mean - shear_field_grid_true, shear_field_grid - shear_field_grid_true]

        plt_data_m = [shear_field_grid_err, num_in_grid]

        plt_data_coms.append(plt_data_com)
        plt_data_ms.append(plt_data_m)

    img = Image_Plot(fig_x=5, fig_y=5 / 4 * 3)
    img.subplots(len(grid_nx), 5)

    titles = ["True", "Mean in grid", "Measured", "Mean - True", "Measured - True"]
    for tag in range(len(grid_nx)):

        nx = grid_nx[tag]
        inverse = range(nx - 1, -1, -1)

        grid_min = min([plt_data_coms[tag][i].min() for i in range(3)])
        grid_max = min([plt_data_coms[tag][i].max() for i in range(3)])

        diff_min = min([plt_data_coms[tag][i].min() for i in range(3,5)])
        diff_max = min([plt_data_coms[tag][i].max() for i in range(3,5)])

        for i in range(len(plt_data_coms[tag])):
            if i < 3:
                vmin, vmax = grid_min, grid_max
            else:
                vmin, vmax = diff_min, diff_max
            fig = img.axs[tag][i].imshow(plt_data_coms[tag][i][inverse],vmin=vmin, vmax=vmax, cmap="jet")
            img.figure.colorbar(fig, ax=img.axs[tag][i])
            img.axs[tag][i].set_title(titles[i])
    img.save_img(parent_path + "pic/grid_com.png")
    img.close_img()

    img = Image_Plot(fig_x=5, fig_y=5 / 4 * 3)
    img.subplots(len(grid_nx), 2)

    titles = ["$\sigma$", "Num"]
    for tag in range(len(grid_nx)):
        nx = grid_nx[tag]
        inverse = range(nx - 1, -1, -1)
        for i in range(len(plt_data_ms[tag])):
            fig = img.axs[tag][i].imshow(plt_data_ms[tag][i][inverse],cmap="jet")
            img.figure.colorbar(fig, ax=img.axs[tag][i])
            img.axs[tag][i].set_title(titles[i])

    img.save_img(parent_path + "pic/grid_m.png")
    img.close_img()
else:
    cmd2 = argv[2]

    result_path = parent_path + "result/grid_filter.hdf5"
    # need 12 threads
    # arcmin
    max_radius = [2, 3, 4]
    # sigma = radius* ratio, 100 means no filter
    sigma_ratio = [0.3, 0.5, 1, 2000]

    # grid num
    nx, ny = 60, 60
    pixel_scale = delta_ra/nx

    # arcmin
    ra_bin = numpy.linspace(-delta_ra / 2, delta_ra / 2, nx + 1)
    dec_bin = numpy.linspace(-delta_dec / 2, delta_dec / 2, ny + 1)

    if cmd2 == "cal":
        # 0 ~ 3 radius[0], 3 ~ 6 radius[1]...
        my_radius, my_sigma = divmod(rank, len(sigma_ratio))
        wei_radius = max_radius[my_radius]
        wei_sigma = sigma_ratio[my_sigma]*wei_radius
        print(rank, my_radius, my_sigma)

        # grid center
        grid_ra = numpy.zeros((ny, nx))
        grid_dec = numpy.zeros((ny, nx))

        # measured shear
        shear_field_grid = numpy.zeros((ny, nx))
        shear_field_grid_err = numpy.zeros((ny, nx))

        # mean shear
        shear_field_grid_mean = numpy.zeros((ny, nx))
        shear_field_grid_mean_w = numpy.zeros((ny, nx))

        # galaxy number
        num_in_grid = numpy.zeros((ny, nx))

        if rank == 0:
            h5f = h5py.File(result_path, "w")
            h5f.close()
        comm.Barrier()

        t1 = time.time()
        for i in range(ny):
            cen_dec = (dec_bin[i] + dec_bin[i + 1]) / 2
            for j in range(nx):
                cen_ra = (ra_bin[j] + ra_bin[j + 1]) / 2

                grid_ra[i, j] = cen_ra
                grid_dec[i, j] = cen_dec

                dx = ra - cen_ra
                dy = dec - cen_dec
                pts_radius = numpy.sqrt(dx ** 2 + dy ** 2)

                idx = pts_radius <= wei_radius

                num_in_grid[i, j] = idx.sum()

                # weight
                weight = MCMC_program.gauss_weight(wei_sigma, pts_radius[idx])
                if sigma_ratio[my_sigma] > 10:
                    weight = numpy.ones_like(weight)

                mg_grid = mg1[idx]*weight
                mnu_grid = mnu1[idx]*weight

                gh, gh_sig = fq.fmin_g_new(mg_grid, mnu_grid, 8)[:2]
                shear_field_grid[i, j] = gh
                shear_field_grid_err[i, j] = gh_sig

                shear_field_grid_mean[i, j] = shear_true[idx].mean()
                shear_field_grid_mean_w[i, j] = numpy.sum(shear_true[idx]*weight)/weight.sum()

        shear_field_grid_true = MCMC_program.shear_slope(parameters, grid_ra, grid_dec)

        t2 = time.time()
        comm.Barrier()

        for i in range(cpus):
            if i == rank:
                h5f = h5py.File(result_path)
                if rank == 0:
                    h5f["/ra"] = grid_ra
                    h5f["/dec"] = grid_dec
                    h5f["/shear_true"] = shear_field_grid_true

                h5f["/shear/raidus_%d_sigma_%d"%(my_radius, my_sigma)] = shear_field_grid
                h5f["/shear_err/raidus_%d_sigma_%d"%(my_radius, my_sigma)] = shear_field_grid_err

                h5f["/num/raidus_%d_sigma_%d"%(my_radius, my_sigma)] = num_in_grid

                h5f["/shear_mean/raidus_%d_sigma_%d"%(my_radius, my_sigma)] = shear_field_grid_mean
                h5f["/shear_mean_wight/raidus_%d_sigma_%d"%(my_radius, my_sigma)] = shear_field_grid_mean_w
                h5f.close()
            comm.Barrier()

        print("%.2f sec"%(t2-t1))

    else:
        if rank == 0:
            h5f = h5py.File(result_path, "r")

            shear_true = h5f["/shear_true"].value

            for ir in range(len(max_radius)):

                pic_nm_shear = parent_path + "pic/grid_filter_shear_%d.png"%ir
                pic_nm_diff = parent_path + "pic/grid_filter_diff_%d.png"%ir
                pic_nm_num = parent_path + "pic/grid_filter_num_%d.png"%ir

                pic_nm_shear_in = parent_path + "pic/grid_filter_shear_%d_in.png"%ir
                pic_nm_diff_in = parent_path + "pic/grid_filter_diff_%d_in.png"%ir
                pic_nm_num_in = parent_path + "pic/grid_filter_num_%d_in.png"%ir

                dn = int(max_radius[ir]/pixel_scale)

                # plot the true and measured shear field
                img = Image_Plot(fig_x=5, fig_y=5/4*3)
                igy, igx = 3, len(sigma_ratio)
                img.subplots(igy, igx)

                titles = ["True", "Mean", "Measured"]

                plt_datas = []
                vmax, vmin = 0, 0
                for ix in range(igx):

                    plt_data = [shear_true,
                                h5f["/shear_mean/raidus_%d_sigma_%d"%(ir, ix)].value,
                                h5f["/shear/raidus_%d_sigma_%d"%(ir, ix)].value]
                    plt_datas.append(plt_data)

                    sub_vmax = max([plt_data[i].max() for i in range(igy)])
                    if vmax < sub_vmax:
                        vmax = sub_vmax

                    sub_vmin = min([plt_data[i].min() for i in range(igy)])
                    if vmin > sub_vmin:
                        vmin = sub_vmin

                for ix in range(igx):
                    for iy in range(igy):
                        fig = img.axs[iy][ix].imshow(plt_datas[ix][iy], vmin=vmin, vmax=vmax, cmap="jet")
                        img.figure.colorbar(fig, ax=img.axs[iy][ix])
                        img.axs[iy][ix].set_title(titles[iy])
                img.save_img(pic_nm_shear)
                img.close_img()
                # plot the true and measured shear field ZOOM IN
                img = Image_Plot(fig_x=5, fig_y=5 / 4 * 3)
                igy, igx = 3, len(sigma_ratio)
                img.subplots(igy, igx)

                titles = ["True", "Mean", "Measured"]

                plt_datas = []
                vmax, vmin = 0, 0
                for ix in range(igx):

                    plt_data = [shear_true[dn:nx-dn, dn:nx-dn],
                                h5f["/shear_mean/raidus_%d_sigma_%d" % (ir, ix)].value[dn:nx-dn, dn:nx-dn],
                                h5f["/shear/raidus_%d_sigma_%d" % (ir, ix)].value[dn:nx-dn, dn:nx-dn]]
                    plt_datas.append(plt_data)

                    sub_vmax = max([plt_data[i].max() for i in range(igy)])
                    if vmax < sub_vmax:
                        vmax = sub_vmax

                    sub_vmin = min([plt_data[i].min() for i in range(igy)])
                    if vmin > sub_vmin:
                        vmin = sub_vmin

                for ix in range(igx):
                    for iy in range(igy):
                        fig = img.axs[iy][ix].imshow(plt_datas[ix][iy], vmin=vmin, vmax=vmax, cmap="jet")
                        img.figure.colorbar(fig, ax=img.axs[iy][ix])
                        img.axs[iy][ix].set_title(titles[iy])
                img.save_img(pic_nm_shear_in)
                img.close_img()


#####################################################################################################
                # plot the difference
                img = Image_Plot(fig_x=5, fig_y=5/4*3)
                igy, igx = 2, len(sigma_ratio)
                img.subplots(igy, igx)

                titles = ["Mean - True", "Measured - True"]

                plt_datas = []
                vmax_1, vmin_1 = 0, 0
                vmax_2, vmin_2 = 0, 0
                for ix in range(igx):

                    plt_data = [h5f["/shear_mean/raidus_%d_sigma_%d"%(ir, ix)].value - shear_true,
                                h5f["/shear/raidus_%d_sigma_%d"%(ir, ix)].value - shear_true]
                    plt_datas.append(plt_data)

                    sub_vmin_1 = plt_data[0].min()
                    sub_vmin_2 = plt_data[1].min()

                    sub_vmax_1 = plt_data[0].max()
                    sub_vmax_2 = plt_data[1].max()

                    if vmax_1 < sub_vmax_1:
                        vmax_1 = sub_vmax_1
                    if vmax_2 < sub_vmax_2:
                        vmax_2 = sub_vmax_2

                    if vmin_1 > sub_vmin_1:
                        vmin_1 = sub_vmin_1
                    if vmin_2 > sub_vmin_2:
                        vmin_2 = sub_vmin_2

                for ix in range(igx):
                    for iy in range(igy):
                        if iy == 0:
                            vmax, vmin = vmax_1, vmin_1
                        else:
                            vmax, vmin = vmax_2, vmin_2
                        fig = img.axs[iy][ix].imshow(plt_datas[ix][iy],vmin=vmin, vmax=vmax, cmap="jet")
                        img.figure.colorbar(fig, ax=img.axs[iy][ix])
                        img.axs[iy][ix].set_title(titles[iy])
                img.save_img(pic_nm_diff)
                img.close_img()

                # plot the difference ZOON IN
                img = Image_Plot(fig_x=5, fig_y=5 / 4 * 3)
                igy, igx = 2, len(sigma_ratio)
                img.subplots(igy, igx)

                titles = ["Mean - True", "Measured - True"]

                plt_datas = []
                vmax_1, vmin_1 = 0, 0
                vmax_2, vmin_2 = 0, 0
                for ix in range(igx):

                    plt_data = [h5f["/shear_mean/raidus_%d_sigma_%d" % (ir, ix)].value[dn:nx-dn, dn:nx-dn] - shear_true[dn:nx-dn, dn:nx-dn],
                                h5f["/shear/raidus_%d_sigma_%d" % (ir, ix)].value[dn:nx-dn, dn:nx-dn] - shear_true[dn:nx-dn, dn:nx-dn]]
                    plt_datas.append(plt_data)

                    sub_vmin_1 = plt_data[0].min()
                    sub_vmin_2 = plt_data[1].min()

                    sub_vmax_1 = plt_data[0].max()
                    sub_vmax_2 = plt_data[1].max()

                    if vmax_1 < sub_vmax_1:
                        vmax_1 = sub_vmax_1
                    if vmax_2 < sub_vmax_2:
                        vmax_2 = sub_vmax_2

                    if vmin_1 > sub_vmin_1:
                        vmin_1 = sub_vmin_1
                    if vmin_2 > sub_vmin_2:
                        vmin_2 = sub_vmin_2

                for ix in range(igx):
                    for iy in range(igy):
                        if iy == 0:
                            vmax, vmin = vmax_1, vmin_1
                        else:
                            vmax, vmin = vmax_2, vmin_2
                        fig = img.axs[iy][ix].imshow(plt_datas[ix][iy], vmin=vmin, vmax=vmax, cmap="jet")
                        img.figure.colorbar(fig, ax=img.axs[iy][ix])
                        img.axs[iy][ix].set_title(titles[iy])
                img.save_img(pic_nm_diff_in)
                img.close_img()

#################################################################################################
                # plot the Number and sigma
                img = Image_Plot(fig_x=5, fig_y=5/4*3)
                igy, igx = 2, len(sigma_ratio)
                img.subplots(igy, igx)

                titles = ["$\sigma$", "Num"]

                plt_datas = []
                vmax_1, vmin_1 = 0, 0
                vmax_2, vmin_2 = 0, 0
                for ix in range(igx):

                    plt_data = [h5f["/shear_err/raidus_%d_sigma_%d"%(ir, ix)].value,
                                h5f["/num/raidus_%d_sigma_%d"%(ir, ix)].value]
                    plt_datas.append(plt_data)

                    sub_vmin_1 = plt_data[0].min()
                    sub_vmin_2 = plt_data[1].min()

                    sub_vmax_1 = plt_data[0].max()
                    sub_vmax_2 = plt_data[1].max()

                    if vmax_1 < sub_vmax_1:
                        vmax_1 = sub_vmax_1
                    if vmax_2 < sub_vmax_2:
                        vmax_2 = sub_vmax_2

                    if vmin_1 > sub_vmin_1:
                        vmin_1 = sub_vmin_1
                    if vmin_2 > sub_vmin_2:
                        vmin_2 = sub_vmin_2

                for ix in range(igx):
                    for iy in range(igy):
                        if iy == 0:
                            vmax, vmin = vmax_1, vmin_1
                        else:
                            vmax, vmin = vmax_2, vmin_2
                        fig = img.axs[iy][ix].imshow(plt_datas[ix][iy], vmin=vmin, vmax=vmax, cmap="jet")
                        img.figure.colorbar(fig, ax=img.axs[iy][ix])
                        img.axs[iy][ix].set_title(titles[iy])
                img.save_img(pic_nm_num)
                img.close_img()

                # plot the Number and sigma ZOON IN
                img = Image_Plot(fig_x=5, fig_y=5 / 4 * 3)
                igy, igx = 2, len(sigma_ratio)
                img.subplots(igy, igx)

                titles = ["$\sigma$", "Num"]

                plt_datas = []
                vmax_1, vmin_1 = 0, 0
                vmax_2, vmin_2 = 0, 0
                for ix in range(igx):

                    plt_data = [h5f["/shear_err/raidus_%d_sigma_%d" % (ir, ix)].value[dn:nx-dn, dn:nx-dn],
                                h5f["/num/raidus_%d_sigma_%d" % (ir, ix)].value[dn:nx-dn, dn:nx-dn]]
                    plt_datas.append(plt_data)

                    sub_vmin_1 = plt_data[0].min()
                    sub_vmin_2 = plt_data[1].min()

                    sub_vmax_1 = plt_data[0].max()
                    sub_vmax_2 = plt_data[1].max()

                    if vmax_1 < sub_vmax_1:
                        vmax_1 = sub_vmax_1
                    if vmax_2 < sub_vmax_2:
                        vmax_2 = sub_vmax_2

                    if vmin_1 > sub_vmin_1:
                        vmin_1 = sub_vmin_1
                    if vmin_2 > sub_vmin_2:
                        vmin_2 = sub_vmin_2

                for ix in range(igx):
                    for iy in range(igy):
                        if iy == 0:
                            vmax, vmin = vmax_1, vmin_1
                        else:
                            vmax, vmin = vmax_2, vmin_2
                        fig = img.axs[iy][ix].imshow(plt_datas[ix][iy], vmin=vmin, vmax=vmax, cmap="jet")
                        img.figure.colorbar(fig, ax=img.axs[iy][ix])
                        img.axs[iy][ix].set_title(titles[iy])
                img.save_img(pic_nm_num_in)
                img.close_img()

            h5f.close()











