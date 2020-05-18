import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
from mpi4py import MPI
import h5py
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import component_fit


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


fq = Fourier_Quad(12, 1234)

data_path = argv[1]
xy_bin_num, radius_bin_num = int(argv[2]), int(argv[3])
shear_scale = float(argv[4])


pic_path = data_path + "/multipole_pic_%.1f"%shear_scale

if rank == 0:
    if not os.path.exists(pic_path):
        os.makedirs(pic_path)


scale = 1000

data_nm = [["noisy_cpp"]]
# data_nm = [["noise_free"], ["noisy_cpp"], ["cross_term"], ["noise_residual"], ["cross_term", "noise_residual"]]

h5f = h5py.File(data_path + "/shear.hdf5", "r")
g1t = h5f["/g1"][()]
g2t = h5f["/g2"][()]
h5f.close()



for tag, nm in enumerate(data_nm):
    pic_nm = "-".join(nm)
    for sub_tag, sub_nm in enumerate(nm):

        if rank == 0:
            print(pic_nm)
        if sub_tag == 0:
            h5f = h5py.File(data_path + "/data_%s_epsf_%d.hdf5"%(sub_nm, rank), "r")
            data = h5f["/data"][()]/scale
            mg1 = data[:,0]
            mg2 = data[:,1]
            mn = data[:,2]
            mu = data[:,3]
            h5f.close()
        else:
            h5f = h5py.File(data_path + "/data_%s_epsf_%d.hdf5"%(sub_nm, rank), "r")
            data = h5f["/data"][()]/scale
            mg1 = mg1 + data[:,0]
            mg2 = mg2 + data[:,1]
            mn = mn + data[:,2]
            mu = mu + data[:,3]
            h5f.close()


    # num, xgrid, ygrid, radius_grid = component_fit.get_bingrid(mg1, mg2, xy_bin_num, 1, 0.3, 99.7)[:4]
    #
    # dpl, radius_bin, radius_mask, mean_of_annuli = component_fit.get_dipole(num, radius_grid, radius_bin_num)
    #
    # qpl, dpl_fit, sin_theta, cos_theta = component_fit.get_quadrupole(dpl, xgrid, ygrid, radius_bin, radius_bin_num)
    #
    # qpl_fit, sin_2theta, cos_2theta = component_fit.fit_quadrupole(qpl,xgrid, ygrid, radius_bin, radius_bin_num)

    # chisq1 = fq.get_chisq_range(mg1, mnu1, 10, gh)[1]
    # chisq2 = fq.get_chisq_range(mg2, mnu2, 10, gh)[1]

    mnu1 = mn + mu
    mnu2 = mn - mu

    pixel_num = 31
    dpl_stack = numpy.zeros((pixel_num*xy_bin_num, xy_bin_num*pixel_num))
    dpl_fit_stack = numpy.zeros((pixel_num*xy_bin_num, xy_bin_num*pixel_num))
    qpl_stack = numpy.zeros((pixel_num*xy_bin_num, xy_bin_num*pixel_num))
    qpl_fit_stack = numpy.zeros((pixel_num*xy_bin_num, xy_bin_num*pixel_num))

    gh1 = numpy.linspace(g1t[rank]-0.05, g1t[rank]+0.05, pixel_num)
    gh2 = numpy.linspace(g2t[rank]-0.05, g2t[rank]+0.05, pixel_num)

    for i in range(pixel_num):
        mg1_sym = mg1 - gh1[i] * mnu1 * shear_scale

        for j in range(pixel_num):
            mg2_sym = mg2 - gh2[j] * mnu2 * shear_scale

            num_sym, xgrid_sym, ygrid_sym, radius_grid_sym = component_fit.get_bingrid(mg1_sym, mg2_sym, xy_bin_num, 1, 0.3, 99.7)[:4]

            dpl_sym, radius_bin_sym, radius_mask_sym, mean_of_annuli_sym = component_fit.get_dipole(num_sym, radius_grid_sym, radius_bin_num)

            qpl_sym, dpl_sym_fit, sin_theta_sym, cos_theta_sym = component_fit.get_quadrupole(dpl_sym, xgrid_sym, ygrid_sym, radius_bin_sym, radius_bin_num)

            qpl_sym_fit, sin_2theta_sym, cos_2theta_sym = component_fit.fit_quadrupole(qpl_sym,xgrid_sym, ygrid_sym, radius_bin_sym, radius_bin_num)

            dpl_stack[i*xy_bin_num: (i+1)*xy_bin_num, j*xy_bin_num:(j+1)*xy_bin_num] = dpl_sym
            dpl_fit_stack[i*xy_bin_num: (i+1)*xy_bin_num, j*xy_bin_num:(j+1)*xy_bin_num] = dpl_sym_fit

            qpl_stack[i*xy_bin_num: (i+1)*xy_bin_num, j*xy_bin_num:(j+1)*xy_bin_num] = qpl_sym
            qpl_fit_stack[i*xy_bin_num: (i+1)*xy_bin_num, j*xy_bin_num:(j+1)*xy_bin_num] = qpl_sym_fit


    img = Image_Plot(fig_x=15, fig_y=15)
    img.subplots(1, 1)
    fig = img.axs[0][0].imshow(dpl_stack)
    img.figure.colorbar(fig, ax=img.axs[0][0])
    img.del_ticks(0, 0, [0, 1])
    pic_name = pic_path + "/%s_%d_compare_2d_dpl.png" % (pic_nm, rank)
    img.save_img(pic_name)
    img.close_img()


    img = Image_Plot(fig_x=15, fig_y=15)
    img.subplots(1, 1)
    fig = img.axs[0][0].imshow(dpl_fit_stack)
    img.figure.colorbar(fig, ax=img.axs[0][0])
    img.del_ticks(0, 0, [0, 1])
    pic_name = pic_path + "/%s_%d_compare_2d_dpl_fit.png" % (pic_nm, rank)
    img.save_img(pic_name)
    img.close_img()


    img = Image_Plot(fig_x=15, fig_y=15)
    img.subplots(1, 1)
    fig = img.axs[0][0].imshow(qpl_stack)
    img.figure.colorbar(fig, ax=img.axs[0][0])
    img.del_ticks(0, 0, [0, 1])
    pic_name = pic_path + "/%s_%d_compare_2d_qpl.png" % (pic_nm, rank)
    img.save_img(pic_name)
    img.close_img()


    img = Image_Plot(fig_x=15, fig_y=15)
    img.subplots(1, 1)
    fig = img.axs[0][0].imshow(qpl_fit_stack)
    img.figure.colorbar(fig, ax=img.axs[0][0])
    img.del_ticks(0, 0, [0, 1])
    pic_name = pic_path + "/%s_%d_compare_2d_qpl_fit.png" % (pic_nm, rank)
    img.save_img(pic_name)
    img.close_img()

    # for i in range(4):
    #     mpl = mpl_stack[i*xy_bin_num:(i+1)*xy_bin_num]
    #     fig = img.axs[i][0].imshow(mpl)
    #     img.figure.colorbar(fig, ax=img.axs[i][1])
    #     img.del_ticks(i, 0, [0, 1])
    #
    # pic_name = pic_path + "/%s_%d_compare_vary_g.png" % (pic_nm, rank)
    # img.save_img(pic_name)
    # img.close_img()

