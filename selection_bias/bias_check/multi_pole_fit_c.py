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


fq = Fourier_Quad(12,1234)

data_path = argv[1]
xy_bin_num, radius_bin_num = int(argv[2]), int(argv[3])
shear_scale = float(argv[4])


pic_path = data_path + "/multipole_pic_%.1f"%shear_scale

if rank == 0:
    if not os.path.exists(pic_path):
        os.makedirs(pic_path)


scale = 1000

data_nm = [["noise_free"], ["noisy_cpp"], ["cross_term"], ["noise_residual"], ["cross_term", "noise_residual"]]

h5f = h5py.File(data_path + "/shear.hdf5", "r")
g1t = h5f["/g1"][()]
g2t = h5f["/g2"][()]
h5f.close()

gh = numpy.linspace(-0.1, 0.1, 21)

for tag, nm in enumerate(data_nm):
    pic_nm = "-".join(nm)
    for sub_tag, sub_nm in enumerate(nm):

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


    num, xgrid, ygrid, radius_grid = component_fit.get_bingrid(mg1, mg2, xy_bin_num, 1, 0.3, 99.7)[:4]

    dpl, radius_bin, radius_mask, mean_of_annuli = component_fit.get_dipole(num, radius_grid, radius_bin_num)

    qpl, dpl_fit, sin_theta, cos_theta = component_fit.get_quadrupole(dpl, xgrid, ygrid, radius_bin, radius_bin_num)

    qpl_fit, sin_2theta, cos_2theta = component_fit.fit_quadrupole(qpl,xgrid, ygrid, radius_bin, radius_bin_num)


    mnu1 = mn + mu
    mnu2 = mn - mu

    mg1_sym = mg1 - g1t[rank]*mnu1*shear_scale
    mg2_sym = mg2 - g2t[rank]*mnu2*shear_scale

    num_sym, xgrid_sym, ygrid_sym, radius_grid_sym = component_fit.get_bingrid(mg1_sym, mg2_sym, xy_bin_num, 1, 0.3, 99.7)[:4]

    dpl_sym, radius_bin_sym, radius_mask_sym, mean_of_annuli_sym = component_fit.get_dipole(num_sym, radius_grid_sym, radius_bin_num)

    qpl_sym, dpl_fit_sym, sin_theta_sym, cos_theta_sym = component_fit.get_quadrupole(dpl_sym, xgrid_sym, ygrid_sym, radius_bin_sym, radius_bin_num)

    qpl_sym_fit, sin_2theta_sym, cos_2theta_sym = component_fit.fit_quadrupole(qpl_sym,xgrid_sym, ygrid_sym, radius_bin_sym, radius_bin_num)


    chisq1 = fq.get_chisq_range(mg1, mnu1, 10, gh)[1]
    chisq2 = fq.get_chisq_range(mg2, mnu2, 10, gh)[1]

    chisq1_sym = fq.get_chisq_range(mg1_sym, mnu1, 10, gh)[1]
    chisq2_sym = fq.get_chisq_range(mg2_sym, mnu2, 10, gh)[1]

    numpy.savez(pic_path + "/cache_%s_%d.npz" % (pic_nm, rank), num, dpl, qpl, dpl_fit, qpl_fit)
    numpy.savez(pic_path + "/cache_%s_sym_%d.npz" % (pic_nm, rank), num_sym, dpl_sym, qpl_sym, dpl_fit_sym, qpl_sym_fit)
    numpy.savez(pic_path + "/cache_%s_chisq_%d.npz" % (pic_nm, rank), chisq1, chisq2, chisq1_sym, chisq2_sym)

    # #################################################################################
    # # hist of data
    #
    # img = Image_Plot(fig_x=6, fig_y=5)
    # img.subplots(2, 3)
    #
    # plot_data = [[num, dpl, dpl_fit], [0,qpl, qpl_fit]]
    # titles = [["G1-G2-hist", "dipole", "dipole-fit"], ["$\chi^2$", "quadrupole","quadrupole-fit"]]
    #
    # img.axs[1][0].plot(gh, chisq1, label="$\chi^2_{g1}$,g1=%.3f" % g1t[rank])
    # img.axs[1][0].plot(gh, chisq2, label="$\chi^2_{g2}$,g2=%.3f" % g2t[rank])
    # img.set_label(1, 0, 0, "$\chi^2$")
    # img.set_label(1, 0, 1, "$\hat{g}$")
    # img.axs[1][0].legend()
    #
    # for i in range(2):
    #     if i == 0:
    #         st = 0
    #     else:
    #         st = 1
    #     for j in range(st,3):
    #         fig = img.axs[i][j].imshow(plot_data[i][j])
    #         img.figure.colorbar(fig, ax=img.axs[i][j])
    #         img.del_ticks(i,j,[0,1])
    #         img.set_label(i,j,0, "+  G1  -")
    #         img.set_label(i,j,1, "-  G2  +")
    #
    #     for j in range(3):
    #         img.axs[i][j].set_title(titles[i][j])
    #
    # pic_name = pic_path + "/%s_%d.png"%(pic_nm, rank)
    # img.save_img(pic_name)
    # img.close_img()
    #
    #
    #
    # #################################################################################
    # # hist of PDF_SYM_data
    #
    # img = Image_Plot(fig_x=6, fig_y=5)
    # img.subplots(2, 3)
    #
    # plot_data = [[num_sym, dpl_sym, dpl_fit_sym], [0, qpl_sym, qpl_sym_fit]]
    # titles =[["PDF_SYM-G1-G2-hist", "PDF_SYM-dipole", "PDF_SYM-dipole-fit"],
    #          ["PDF_SYM-$\chi^2$", "PDF_SYM-quadrupole", "PDF_SYM-quadrupole-fit"]]
    #
    # img.axs[1][0].plot(gh, chisq1_sym, label="$\chi^2_{g1}$,g1=%.3f" % g1t[rank])
    # img.axs[1][0].plot(gh, chisq2_sym, label="$\chi^2_{g2}$,g2=%.3f" % g2t[rank])
    #
    # img.set_label(1, 0, 0, "$\chi^2$")
    # img.set_label(1, 0, 1, "$\hat{g}$")
    # img.axs[1][0].legend()
    #
    # for i in range(2):
    #     if i == 0:
    #         st = 0
    #     else:
    #         st = 1
    #     for j in range(st, 3):
    #         fig = img.axs[i][j].imshow(plot_data[i][j])
    #         img.figure.colorbar(fig, ax=img.axs[i][j])
    #         img.del_ticks(i, j, [0, 1])
    #         img.set_label(i, j, 0, "+  G1  -")
    #         img.set_label(i, j, 1, "-  G2  +")
    #
    #     for j in range(3):
    #         img.axs[i][j].set_title(titles[i][j])
    #
    # pic_name = pic_path + "/%s_%d_sym.png" % (pic_nm, rank)
    # img.save_img(pic_name)
    # img.close_img()


    #################################################################################
    # compare

    img = Image_Plot(fig_x=5, fig_y=4,xpad=0.25, ypad=0.25)

    img.subplots(2, 3)


    titles = [["$\chi^2$", "dipole-fit", "quadrupole-fit"],
              ["$\chi^2-SYM$", "dipole-fit_SYM", "quadruple-fit_SYM"]]

    img.axs[0][0].plot(gh, chisq1, label="$\chi^2_{g1}$,g1=%.3f" % g1t[rank])
    img.axs[0][0].plot(gh, chisq2, label="$\chi^2_{g2}$,g2=%.3f" % g2t[rank])

    img.axs[1][0].plot(gh, chisq1_sym, label="$\chi^2_{g1}$,g1=%.3f" % g1t[rank])
    img.axs[1][0].plot(gh, chisq2_sym, label="$\chi^2_{g2}$,g2=%.3f" % g2t[rank])

    img.set_label(0, 0, 0, "$\chi^2$")
    img.set_label(0, 0, 1, "$\hat{g}$")
    img.axs[0][0].legend()

    img.set_label(1, 0, 0, "$\chi^2-SYM$")
    img.set_label(1, 0, 1, "$\hat{g}$")
    img.axs[1][0].legend()

    dpl_fit = numpy.nan_to_num(dpl_fit)
    dpl_fit_sym = numpy.nan_to_num(dpl_fit_sym)

    qpl_fit = numpy.nan_to_num(qpl_fit)
    qpl_sym_fit = numpy.nan_to_num(qpl_sym_fit)

    vmax_dpl = max(numpy.abs(dpl_fit).max(), numpy.abs(dpl_fit_sym).max())
    vmax_qpl = max(numpy.abs(qpl_fit).max(), numpy.abs(qpl_sym_fit).max())

    fig = img.axs[0][1].imshow(dpl_fit, vmin=-vmax_dpl, vmax=vmax_dpl)
    img.figure.colorbar(fig, ax=img.axs[0][1])
    fig = img.axs[1][1].imshow(dpl_fit_sym, vmin=-vmax_dpl, vmax=vmax_dpl)
    img.figure.colorbar(fig, ax=img.axs[1][1])

    fig = img.axs[0][2].imshow(qpl_fit, vmin=-vmax_qpl, vmax=vmax_qpl)
    img.figure.colorbar(fig, ax=img.axs[0][2])
    fig = img.axs[1][2].imshow(qpl_sym_fit, vmin=-vmax_qpl, vmax=vmax_qpl)
    img.figure.colorbar(fig, ax=img.axs[1][2])

    for i in range(2):
        for j in range(3):
            if j > 0:
                img.del_ticks(i, j, [0, 1])
                img.set_label(i, j, 0, "+  G1  -")
                img.set_label(i, j, 1, "-  G2  +")

            img.axs[i][j].set_title(titles[i][j])

    pic_name = pic_path + "/%s_%d_compare.png" % (pic_nm, rank)
    img.save_img(pic_name)
    img.close_img()


    # #################################################################################
    # # x & y grid, raidus ....
    # img = Image_Plot()
    # img.subplots(2, 5)
    #
    # titles = [["x-grid", "y-grid", "radius-grid", "radius-bin", "mean_num_annuli"],
    #           ["PDF_SYM-x-grid", "PDF_SYM-y-grid", "PDF_SYM-radius-grid", "PDF_SYM-radius-bin", "PDF_SYM-mean_num_annuli"]]
    # plot_data = [[xgrid, ygrid, radius_grid, radius_mask, mean_of_annuli],
    #              [xgrid_sym, ygrid_sym, radius_grid_sym, radius_mask_sym, mean_of_annuli_sym]]
    # for i in range(2):
    #     for j in range(5):
    #         fig = img.axs[i][j].imshow(plot_data[i][j])
    #         img.figure.colorbar(fig, ax=img.axs[i][j])
    #         img.del_ticks(i,j,[0,1])
    #         img.axs[i][j].set_title(titles[i][j])
    # pic_name = pic_path + "/%s_%d_check1.png"%(pic_nm, rank)
    # img.save_img(pic_name)
    # img.close_img()
    #
    #
    # #################################################################################
    # # sin_theta .....
    # img = Image_Plot()
    # img.subplots(2, 4)
    #
    # titles = [["sin$\\theta$", "cos$\\theta$", "sin$2\\theta$", "cos$2\\theta$"],
    #           ["PDF_SYM-sin$\\theta$", "PDF_SYM-cos$\\theta$", "PDF_SYM-sin$2\\theta$", "PDF_SYM-cos$2\\theta$"]]
    # plot_data = [[sin_theta, cos_theta, sin_2theta, cos_2theta],
    #              [sin_theta_sym, cos_theta_sym, sin_2theta_sym, cos_2theta_sym]]
    #
    # for i in range(2):
    #     for j in range(4):
    #         fig = img.axs[i][j].imshow(plot_data[i][j])
    #         img.figure.colorbar(fig, ax=img.axs[i][j])
    #         img.del_ticks(i,j,[0,1])
    #         img.axs[i][j].set_title(titles[i][j])
    # pic_name = pic_path + "/%s_%d_check2.png"%(pic_nm, rank)
    # img.save_img(pic_name)
    # img.close_img()
