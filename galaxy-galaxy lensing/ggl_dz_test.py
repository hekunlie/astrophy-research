from sys import path,argv
path.append("/home/hklee/work/mylib")
from hk_plot_tool import Image_Plot
import hk_tool_box
import hk_gglensing_tool
import numpy
import h5py
import galsim
import hk_FQlib
from astropy.cosmology import FlatLambdaCDM
from mpi4py import MPI



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


# cosmology
omega_m0 = 0.31
omega_lam0 = 1 - omega_m0
h = 0.6735
H_0 = 100 * h

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)

# Halo parameters
Mass = 3*10 ** 13.5  # M_sun/h
conc = 6  # concentration
len_z = 0.3  # redshift
halo_position = galsim.PositionD(0, 0)  # arcsec
com_dist_len = cosmos.comoving_distance(len_z).value * h  # Mpc/h
print("Lens plane at z = %.2f, %.5f Mpc/h" % (len_z, com_dist_len))


# lens profile
CF = hk_gglensing_tool.Cosmos_flat(omega_m0, 100*h)
CF.NFW((0,0), Mass, conc, len_z)


data_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_6"
test_num = numprocs
task_list = hk_tool_box.alloc([i for i in range(test_num)], numprocs)[rank]
ds_shape = (9, test_num)

sigma_z = float(argv[1])

for dz in [0, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40]:

    ds_result = numpy.zeros(ds_shape)
    rng = numpy.random.RandomState(12312 + rank*100)

    for i in task_list:
        h5f = h5py.File(data_path + "/data/segment_sheared_noisy_cpp.hdf5", "r")
        src_z_ori = h5f["/%d/z" % i][()]
        ra_ori = h5f["/%d/ra" % i][()]
        dec_ori = h5f["/%d/dec" % i][()]
        sep_radius_ori = h5f["/%d/sep_radius" % i][()]
        sep_arcsec_ori = h5f["/%d/sep_radian" % i][()] / numpy.pi * 180 * 3600
        g1_ori = h5f["/%d/g1" % i][()]
        g2_ori = h5f["/%d/g2" % i][()]
        mgt_ori = h5f["/%d/mgt" % i][()]
        mnu1_ori = h5f["/%d/mnu1" % i][()]
        h5f.close()

        idx_ori = src_z_ori > len_z + dz

        src_z_err = rng.normal(0, (1 + src_z_ori) * sigma_z)
        src_z = src_z_ori + src_z_err

        idx = src_z > len_z + dz
        mix_num = idx.sum()

        idx_true = idx_ori & idx
        true_num = idx_true.sum()
        ratio = true_num / mix_num
        print(true_num, mix_num, true_num / mix_num)

        src_z = src_z[idx]
        ra, dec = ra_ori[idx], dec_ori[idx]

        sep_radius = sep_radius_ori[idx]
        sep_arcsec = sep_arcsec_ori[idx]
        g1, g2 = g1_ori[idx], g2_ori[idx]

        mgt = mgt_ori[idx]
        mnu1 = mnu1_ori[idx]

        sigma_crit = CF.get_sigma_crit(src_z)
        count_weight = 1. / (sigma_crit)

        # the true signal
        gt_true = numpy.sqrt(g1_ori[idx_true] ** 2 + g2_ori[idx_true] ** 2)
        sigma_crit_true = CF.get_sigma_crit(src_z_ori[idx_true])
        ds_true = gt_true * sigma_crit_true

        ds_result[0, i] = ratio

        # the pure src signal
        ds_result[1, i] = sep_radius_ori[idx_true].mean()
        ds_result[2, i] = CF.get_shear(sep_arcsec_ori[idx_true].mean(), 0, src_z_ori[idx_true].mean())[2]

        # the true signal measured at the mixed radius
        ds_result[3, i] = sep_radius.mean()
        ds_result[4, i] = CF.get_shear(sep_arcsec.mean(), 0, src_z_ori[idx_true].mean())[2]



        #     img = Image_Plot()
        #     img.subplots(1,2)

        gh, gh_sig = hk_FQlib.find_shear_cpp(mgt, mnu1 / sigma_crit, numpy.ones_like(mgt), bin_num=10, left=-200, right=200,
                                             fit_scale=1)[:2]  # ,fig_ax=img.axs[0][0])[:2]
        ghw, ghw_sig = hk_FQlib.find_shear_cpp(mgt, mnu1 / sigma_crit, count_weight, 10, left=-200, right=200, fit_scale=1)[:2]  # , fig_ax=img.axs[0][1])[:2]

        #     img.show_img()
        ds_result[5, i] = gh
        ds_result[6, i] = gh_sig
        ds_result[7, i] = ghw
        ds_result[8, i] = ghw_sig


        print("%.4f->%.4f.  %.4f(%.4f).  %.4f(%.4f). relative err: %.4f" % (
        ds_true.mean(), ds_true.mean() * ratio, gh, gh_sig, ghw, ghw_sig, (gh_sig - ghw_sig) / ghw_sig))



    # rank 0 will collect all data, the others send data to it.
    if rank > 0:
        # !!!! remember the data type, MPI.DOUBLE, MPI.FLOAT, ...
        # or it will raise an error, Keyerror
        comm.Send([ds_result, MPI.DOUBLE], dest=0, tag=rank)
    else:

        # receive the data from other CPUs
        # !!!! the start points is 1 in range() not 0
        for procs in range(1, numprocs):
            # prepare a buffer for the data, the shape must be the same
            # with that of what the other CPUs send, you have collected them in 'data_sps'
            recvs = numpy.empty(ds_shape, dtype=numpy.double)
            # receive it using the buffer,
            comm.Recv(recvs, source=procs, tag=procs)
            ds_result += recvs

        img = Image_Plot(xpad=0.3, ypad=0.25)
        img.subplots(1, 2)
        img.axs[0][0].plot(ds_result[1], ds_result[2], label="Pure src signal")
        img.axs[0][0].plot(ds_result[3], ds_result[4], label="True signal at the mixed radius")

        img.axs[0][0].errorbar(ds_result[3], ds_result[5], ds_result[6], capsize=3, label="ori")
        img.axs[0][0].errorbar(ds_result[3], ds_result[7], ds_result[8], capsize=3, label="weight")

        img.axs[0][1].errorbar(ds_result[3], ds_result[5] - ds_result[4], ds_result[6], capsize=3, label="ori")
        img.axs[0][1].errorbar(ds_result[3], ds_result[7] - ds_result[4], ds_result[8], capsize=3, label="weight")

        img.axs[0][1].plot(ds_result[3], ds_result[4]*ds_result[0], label="Should Be")

        img.axs[0][0].set_yscale("log")
        for i in range(2):
            img.axs[0][i].set_xscale("log")
            img.axs[0][i].legend()
        img.set_label(0,0,0,"$\Delta\Sigma$")
        img.set_label(0,1,0,"$\delta\Delta\Sigma$")
        img.set_label(0,0,1,"R [Mpc/h]")
        img.set_label(0,1,1,"R [Mpc/h]")
        # img.show_img()
        img.save_img("./test_sigmaz_%.2f_dz_%.2f.png"%(sigma_z, dz))

        numpy.savez("./test_sigmaz_%.2f_dz_%.2f.npz"%(sigma_z, dz), ds_result)
    comm.Barrier()
