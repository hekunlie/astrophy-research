import os
from sys import path
path.append("/home/hklee/work/mylib")
from hk_plot_tool import Image_Plot
import hk_tool_box
import hk_gglensing_tool
import numpy
import h5py
import galsim
import hk_FQlib
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
import time
from mpi4py import MPI

#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

data_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z"


h5f = h5py.File(data_path + "/params/stack_sheared_para.hdf5","r")
src_ra = h5f["/ra"][()]
src_dec = h5f["/dec"][()]
src_z = h5f["/z"][()]
h5f.close()

h5f = h5py.File(data_path + "/data/stack_sheared_data_noise_free.hdf5","r")
shear_est_nf = h5f["/data"][()]
h5f.close()
h5f = h5py.File(data_path + "/data/stack_sheared_data_noisy_cpp.hdf5","r")
shear_est_ny = h5f["/data"][()]
h5f.close()


# cosmology
omega_m0 = 0.31
omega_lam0 = 1 - omega_m0
h = 0.6735
C_0_hat = 2.99792458
H_0 = 100 * h
coeff = 1000 * C_0_hat / h

coeff_crit = C_0_hat ** 2 / 4 / numpy.pi / 6.674

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)

# Halo parameters
Mass = 3*10 ** 13.5  # M_sun/h
conc = 6  # concentration
len_z = 0.3  # redshift
halo_position = galsim.PositionD(0, 0)  # arcsec
com_dist_len = cosmos.comoving_distance(len_z).value * h  # Mpc/h

com_dist_src = cosmos.comoving_distance(src_z).value * h  # Mpc/h

nfw_model = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)

crit_sd = 1662895.2081868195*com_dist_src/com_dist_len/(com_dist_src-com_dist_len)/(1+len_z)

# position and separation angle
pos_len = SkyCoord(ra=0*units.deg, dec=0*units.deg,frame="fk5")
pos_src = SkyCoord(ra=src_ra*units.deg, dec=src_dec*units.deg,frame="fk5")

separation_radian = pos_len.separation(pos_src).radian
separation_radius = separation_radian*com_dist_len
print(separation_radius.min(),separation_radius.max())

position_angle = pos_len.position_angle(pos_src).radian

sin_2theta = numpy.sin(2*position_angle)
cos_2theta = numpy.cos(2*position_angle)

sin_4theta = numpy.sin(4*position_angle)
cos_4theta = numpy.cos(4*position_angle)

mgt_nf = shear_est_nf[:,0]*cos_2theta - shear_est_nf[:,1]*sin_2theta
mgx_nf = shear_est_nf[:,0]*sin_2theta + shear_est_nf[:,1]*cos_2theta
mu_nf = shear_est_nf[:,3]*cos_4theta - shear_est_nf[:,4]*sin_4theta
mn_nf = shear_est_nf[:,2]

mgt_ny = shear_est_ny[:,0]*cos_2theta - shear_est_ny[:,1]*sin_2theta
mgx_ny = shear_est_ny[:,0]*sin_2theta + shear_est_ny[:,1]*cos_2theta
mu_ny = shear_est_ny[:,3]*cos_4theta - shear_est_ny[:,4]*sin_4theta
mn_ny = shear_est_ny[:,2]




radius_bin_num = numprocs
radius_bin = hk_tool_box.set_bin_log(0.2, 19.5, radius_bin_num + 1)
radius_bin_tag = [i for i in range(radius_bin_num)]
my_radius_bin_tag = hk_tool_box.alloc(radius_bin_tag, numprocs)[rank]


itemsize = MPI.DOUBLE.Get_size()

if rank == 0:
    # bytes for 10 double elements
    nbytes1 = 12*radius_bin_num*itemsize
    nbytes2 = 12*radius_bin_num*itemsize
    nbytes3 = 3*radius_bin_num*itemsize
else:
    nbytes1 = 0
    nbytes2 = 0
    nbytes3 = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes1, itemsize, comm=comm)
win2 = MPI.Win.Allocate_shared(nbytes2, itemsize, comm=comm)
win3 = MPI.Win.Allocate_shared(nbytes3, itemsize, comm=comm)
# create a numpy array whose data points to the shared block
# buf is the block's address in the memory
buf1, itemsize = win1.Shared_query(0)
buf2, itemsize = win2.Shared_query(0)
buf3, itemsize = win3.Shared_query(0)

Ds_nf = numpy.ndarray(buffer=buf1, dtype='d', shape=(12, radius_bin_num)) # array filled with zero
Ds_ny = numpy.ndarray(buffer=buf2, dtype='d', shape=(12, radius_bin_num))
inform = numpy.ndarray(buffer=buf3, dtype='d', shape=(3, radius_bin_num))

if rank == 0:
    for i in range(radius_bin_num):
        idx1 = separation_radius >= radius_bin[i]
        idx2 = separation_radius < radius_bin[i + 1]
        idx = idx1 & idx2

        inform[0, i] = separation_radius[idx].mean()
        inform[1, i] = separation_radian[idx].mean() / numpy.pi * 180 * 3600
        inform[2, i] = idx.sum()
comm.Barrier()

#  Model
Ds_true = hk_gglensing_tool.get_delta_sigma(nfw_model, com_dist_len, len_z, com_dist_src[0], src_z[0], inform[1])
if rank == 0:
    print(inform[2])

t1 = time.time()

pdf_bin_num = [2, 10, 20]

for i in my_radius_bin_tag:
    idx1 = separation_radius >= radius_bin[i]
    idx2 = separation_radius < radius_bin[i + 1]
    idx = idx1 & idx2

    temp_mgt_nf = mgt_nf[idx] * crit_sd[idx]
    temp_mgx_nf = mgx_nf[idx] * crit_sd[idx]
    temp_mgnu1_nf = (mn_nf[idx] + mu_nf[idx])#/crit_sd[idx]
    temp_mgnu2_nf = (mn_nf[idx] - mu_nf[idx])#/crit_sd[idx]

    for j in range(3):
        result_t = hk_FQlib.find_shear_cpp(temp_mgt_nf, temp_mgnu1_nf, bin_num=pdf_bin_num[j], left=-200, right=200,
                                           chi_gap=20,max_iters=60)[:2]
        result_x = hk_FQlib.find_shear_cpp(temp_mgx_nf, temp_mgnu2_nf, bin_num=pdf_bin_num[j], left=-200, right=200,
                                           chi_gap=20,max_iters=60)[:2]
        st, ed = int(j * 4), int((j + 1) * 4)
        Ds_nf[st:ed, i] = result_t[0], result_t[1], result_x[0], result_x[1]

    temp_mgt_ny = mgt_ny[idx] * crit_sd[idx]
    temp_mgx_ny = mgx_ny[idx] * crit_sd[idx]
    temp_mgnu1_ny = (mn_ny[idx] + mu_ny[idx])#/crit_sd[idx]
    temp_mgnu2_ny = (mn_ny[idx] - mu_ny[idx])#/crit_sd[idx]

    for j in range(3):
        result_t = hk_FQlib.find_shear_cpp(temp_mgt_ny, temp_mgnu1_ny, bin_num=pdf_bin_num[j], left=-200, right=200,
                                           chi_gap=20,max_iters=60)[:2]
        result_x = hk_FQlib.find_shear_cpp(temp_mgx_ny, temp_mgnu2_ny, bin_num=pdf_bin_num[j], left=-200, right=200,
                                           chi_gap=20,max_iters=60)[:2]
        st, ed = int(j * 4), int((j + 1) * 4)
        Ds_ny[st:ed, i] = result_t[0], result_t[1], result_x[0], result_x[1]
t2 = time.time()
comm.Barrier()
print(rank, t2-t1)

if rank == 0:


    result_path = "."
    numpy.savez(result_path+"/result_%d.npz"%numprocs, inform,Ds_nf,Ds_ny)

    img = Image_Plot(xpad=0.25, ypad=0.24)
    img.subplots(len(pdf_bin_num), 4)

    for j in range(len(pdf_bin_num)):
        tag = int(4 * j)
        img.axs[j][0].errorbar(inform[0], Ds_nf[tag], Ds_nf[tag + 1], capsize=3, marker="s", fmt=" ",
                               label="Noise free")
        img.axs[j][0].errorbar(inform[0], Ds_ny[tag], Ds_ny[tag + 1], capsize=3, marker="o", fmt=" ",
                               label="Noisy")

        img.axs[j][1].errorbar(inform[0], Ds_nf[tag + 2], Ds_nf[tag + 3], capsize=3, marker="s", fmt=" ",
                               label="Noise free")
        img.axs[j][1].errorbar(inform[0], Ds_ny[tag + 2], Ds_ny[tag + 3], capsize=3, marker="o", fmt=" ",
                               label="Noisy")

        img.axs[j][0].plot(inform[0], Ds_true, label="Model")

        img.axs[j][2].plot(inform[0], (Ds_nf[tag] - Ds_true),marker="s", label="Noise free")
        img.axs[j][2].plot(inform[0], (Ds_ny[tag] - Ds_true), marker="o",label="Noisy")
        # img.axs[0][2].set_yscale("log")

        img.axs[j][3].plot(inform[0], (Ds_nf[ tag] - Ds_true) / Ds_nf[tag],marker="s", label="Noise free")
        img.axs[j][3].plot(inform[0], (Ds_ny[tag] - Ds_true) / Ds_ny[tag], marker="o",label="Noisy")

        for i in range(4):
            img.set_label(j, i, 1, "Radius Mpc/h")
            img.axs[j][i].set_xscale("log")
            img.axs[j][i].legend()
        img.set_label(j, 0, 0, "$\Delta\Sigma$")
        img.set_label(j, 1, 0, "$\Delta\Sigma_x$")

        img.set_label(j, 2, 0, "$\Delta\Sigma - \Delta\Sigma_{model}$")
        img.set_label(j, 3, 0, "$(\Delta\Sigma - \Delta\Sigma_{model})/Error bar$")

        img.axs[j][0].set_yscale("log")
    img.save_img(result_path + "/signal_comparison_%d.pdf"%numprocs)
    # img.show_img()

    img = Image_Plot(xpad=0.25)
    img.subplots(1, 2)
    for i in range(1, 3):
        tag = int(i * 4)
        img.axs[0][0].plot(inform[0], Ds_nf[tag + 1] / Ds_nf[1],
                           label="Noise free: Error bar %d bins/ %d bins" % (pdf_bin_num[i], pdf_bin_num[0]))
        img.axs[0][1].plot(inform[0], Ds_ny[tag + 1] / Ds_ny[1],
                           label="Noisy: Error bar %d bins/ %d bins" % (pdf_bin_num[i], pdf_bin_num[0]))
    for i in range(2):
        img.axs[0][i].legend()
        img.axs[0][i].set_xscale("log")
        img.set_label(0, i, 0, "Error bar ratio")
        img.set_label(0, i, 1, "Radius Mpc/h")
    img.save_img(result_path + "/err_comparison_%d.pdf"%numprocs)
    # img.show_img()
comm.Barrier()