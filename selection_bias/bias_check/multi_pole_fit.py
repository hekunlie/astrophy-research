import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
from mpi4py import MPI
import h5py
from scipy import optimize
from plot_tool import Image_Plot

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


def dipole_fun(sin_cos, a1, a2):
    return a1 * sin_cos[0] + a2 * sin_cos[1]


def get_bingrid(data_1, data_2, bin_num, mode=0, percen_low=0.5, percen_up=99.5):
    if mode == 0:
        num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [bin_num, bin_num])
    #         print(mode, xbins - ybins)

    elif mode == 1:
        a = max(numpy.abs(numpy.percentile(data_1, percen_low)), numpy.percentile(data_1, percen_up))
        b = max(numpy.abs(numpy.percentile(data_2, percen_low)), numpy.percentile(data_2, percen_up))
        xy_max = max(a, b)
        xy_bin = numpy.linspace(-xy_max, xy_max, bin_num + 1)
        num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [xy_bin, xy_bin])
    #         print(mode, xbins-ybins)
    else:
        xy_max = max(numpy.max(numpy.abs(data_1)), numpy.max(numpy.abs(data_2)))
        xy_bin = numpy.linspace(-xy_max, xy_max, bin_num + 1)
        num, d1_bins, d2_bins = numpy.histogram2d(data_1, data_2, [xy_bin, xy_bin])
    #         print(mode, xbins-ybins)

    xgrid = numpy.zeros((bin_num, bin_num))
    ygrid = numpy.zeros((bin_num, bin_num))

    for i in range(bin_num):
        for j in range(bin_num):
            xgrid[i, j] = (d2_bins[j] + d2_bins[j + 1]) / 2
            ygrid[i, j] = (d1_bins[i] + d1_bins[i + 1]) / 2
    return num, xgrid, ygrid, numpy.sqrt(xgrid ** 2 + ygrid ** 2)


def get_dipole(num, radius, radius_bin_num):
    radius_bin = numpy.linspace(0, radius.max(), radius_bin_num + 1)
    num_dipole = numpy.zeros_like(num)
    raidus_mask = numpy.zeros_like(num)
    for i in range(radius_bin_num):
        idx1 = radius >= radius_bin[i]
        idx2 = radius < radius_bin[i + 1]
        idx = idx1 & idx2
        num_dipole[idx] = num[idx] - num[idx].mean()
        raidus_mask[idx] = i
    return num_dipole, radius_bin, raidus_mask


def get_quadrupole(num_dipole, xgrid, ygrid, radius_bin, radius_bin_num):
    radius_grid = numpy.sqrt(xgrid**2 + ygrid**2)
    sin_theta = ygrid / radius_grid
    cos_theta = xgrid / radius_grid
    num_dipole_fit = numpy.zeros_like(num_dipole)
    num_dipole_fit[:, :] = numpy.nan
    #     print(num_dipole_fit)
    for i in range(radius_bin_num):
        idx1 = radius_grid >= radius_bin[i]
        idx2 = radius_grid < radius_bin[i + 1]
        idx = idx1 & idx2
        if idx.sum() > 5:
            sin_cos = numpy.zeros((2, idx.sum()))
            sin_cos[0] = sin_theta[idx]
            sin_cos[1] = cos_theta[idx]
            res = optimize.curve_fit(dipole_fun, sin_cos, num_dipole[idx])[0]
            #             print(res)
            num_dipole_fit[idx] = dipole_fun(sin_cos, res[0], res[1])
    return num_dipole - num_dipole_fit, num_dipole_fit



data_path = argv[1]
xy_bin_num, radius_bin_num = int(argv[2]), int(argv[3])

pic_path = data_path + "/multipole_pic"

if rank == 0:
    if not os.path.exists(pic_path):
        os.makedirs(pic_path)


scale = 1000

data_nm = [["noise_free"], ["noisy_cpp"], ["cross_term"], ["noise_residual"], ["cross_term", "noise_residual"]]

h5f = h5py.File(data_path + "/shear.hdf5", "r")
g1t = h5f["/g1"][()]
g2t = h5f["/g2"][()]
h5f.close()

for tag, nm in enumerate(data_nm):
    pic_nm = ""
    for sub_tag, sub_nm in enumerate(nm):

        pic_nm += sub_nm
        if rank == 0:
            print(pic_nm)
        if sub_tag == 0:
            h5f = h5py.File(data_path + "/data_%s_%d.hdf5"%(sub_nm, rank), "r")
            data = h5f["/data"][()]/scale
            mg1 = data[:,0]
            mg2 = data[:,1]
            mn = data[:,2]
            mu = data[:,3]
            h5f.close()
        else:
            h5f = h5py.File(data_path + "/data_%s_%d.hdf5"%(sub_nm, rank), "r")
            data = h5f["/data"][()]/scale
            mg1 = mg1 + data[:,0]
            mg2 = mg2 + data[:,1]
            mn = mn + data[:,2]
            mu = mu + data[:,3]
            h5f.close()

    mg1_sym = mg1 - g1t[rank]*(mn + mu)
    mg2_sym = mg2 - g2t[rank]*(mn - mu)

    num, xgrid, ygrid, radius_grid = get_bingrid(mg1, mg2, xy_bin_num, 1, 0.05, 99.95)
    dpl, radius_bin, radius_mask = get_dipole(num, radius_grid, radius_bin_num)
    qpl, dpl_fit = get_quadrupole(dpl, xgrid, ygrid, radius_bin, radius_bin_num)

    num_sym, xgrid_sym, ygrid_sym, radius_grid_sym = get_bingrid(mg1_sym, mg2_sym, xy_bin_num, 1, 0.05, 99.95)
    dpl_sym, radius_bin_sym, radius_mask_sym = get_dipole(num_sym, radius_grid_sym, radius_bin_num)
    qpl_sym, dpl_fit_sym = get_quadrupole(dpl_sym, xgrid_sym, ygrid_sym, radius_bin_sym, radius_bin_num)

    plot_data = [[num, dpl, dpl_fit, qpl],
                 [num_sym, dpl_sym, dpl_fit_sym, qpl_sym]]

    img = Image_Plot()
    img.subplots(2, 4)
    for i in range(2):
        for j in range(4):
            fig = img.axs[i][j].imshow(plot_data[i][j])
            img.figure.colorbar(fig, ax=img.axs[i][j])
            img.del_tick(i,j,[0,1])
    pic_name = pic_path + "/%s_%d.png"%(pic_nm, rank)
    img.save_img(pic_name)
    img.close_img()

    plot_data = [[xgrid, ygrid, radius_grid, radius_mask],
                 [xgrid_sym, ygrid_sym, radius_grid_sym, radius_mask_sym]]

    img = Image_Plot()
    img.subplots(2, 4)
    for i in range(2):
        for j in range(4):
            fig = img.axs[i][j].imshow(plot_data[i][j])
            img.figure.colorbar(fig, ax=img.axs[i][j])
            img.del_tick(i,j,[0,1])
    pic_name = pic_path + "/debug_%s_%d.png"%(pic_nm, rank)
    img.save_img(pic_name)
    img.close_img()