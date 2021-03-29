from sys import path,argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import gglensing_tool
import numpy
import h5py
import galsim
from Fourier_Quad import Fourier_Quad
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
import time
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

mg_bin_num = int(argv[1])

# cosmology
omega_m0 = 0.31
omega_lam0 = 1 - omega_m0
h = 0.6735
C_0_hat = 2.99792458
H_0 = 100*h
coeff = 1000*C_0_hat/h

coeff_crit = C_0_hat**2/4/numpy.pi/6.674

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)


# Halo parameters
Mass = 5*10**12.5 #M_sun/h
conc = 3.5 #concentration
len_z = 0.3 #redshift
halo_position = galsim.PositionD(0,0) #arcsec
com_dist_len = cosmos.comoving_distance(len_z).value*h #Mpc/h
# print("Lens plane at z = %.2f, %.5f Mpc/h"%(len_z, com_dist_len) )
len_pos = SkyCoord(ra=0*units.arcsec, dec=0*units.arcsec,frame="fk5")

# lens profile
nfw = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)

shear_tag = 0

data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_2"

data_path1 = data_path + "/data/sheared_data_%d_noisy_cpp.hdf5" % shear_tag
para_path1 = data_path + "/sheared_para_%d.hdf5" % shear_tag

data_path2 = data_path + "/data/non_sheared_data_%d_noisy_cpp.hdf5" % shear_tag
para_path2 = data_path + "/non_sheared_para_%d.hdf5" % shear_tag

data_path3 = data_path + "/data/non_sheared_data_%d_noisy_cpp.hdf5" % (shear_tag + 5)
para_path3 = data_path + "/non_sheared_para_%d.hdf5" % (shear_tag + 5)


num_s = 1000000

foreground_z_err = numpy.abs(numpy.random.normal(0, 0.2, 1000000)) + 0.31

dilute_case = 5
max_dilute_ratio = 0.4
corr_case = 51
dilute_ratio_list = numpy.linspace(0,max_dilute_ratio,dilute_case)
corr_num = numpy.linspace(0, 500000, corr_case,dtype=numpy.intc)

task_list = [[i,j] for i in range(dilute_case) for j in range(corr_case)]
task_list_sub = tool_box.alloc(task_list, numprocs)[rank]


itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = (corr_case)*(dilute_case+dilute_case)*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)

total_result = numpy.ndarray(buffer=buf1, dtype='d', shape=(dilute_case+dilute_case, corr_case))


comm.Barrier()

for ij in task_list_sub:

    dilute_ratio = dilute_ratio_list[ij[0]]
    num_non = int(num_s * dilute_ratio)
    num_corr = corr_num[ij[1]]

    mg1, mg2, mn, mu, mv, ra, dec, z = gglensing_tool.data_mix(data_path1, para_path1, data_path2, para_path2,
                                                               num_s, num_non, foreground_z_err[:num_non])

    com_dist = cosmos.comoving_distance(z).value * H_0 / 100

    crit_coeff = 1662895.2081868195 * com_dist / com_dist_len / (com_dist - com_dist_len) / (1 + len_z)

    src_pos = SkyCoord(ra=ra * units.arcsec, dec=dec * units.arcsec, frame="fk5")
    position_theta = len_pos.position_angle(src_pos).rad

    cos_2theta = numpy.cos(2 * position_theta)
    sin_2theta = numpy.sin(2 * position_theta)
    cos_4theta = numpy.cos(4 * position_theta)
    sin_4theta = numpy.sin(4 * position_theta)

    mg1r = (mg1 * cos_2theta - mg2 * sin_2theta) * crit_coeff
    mg2r = (mg1 * sin_2theta + mg2 * cos_2theta) * crit_coeff
    mur = mu * cos_4theta - mv * sin_4theta
    mnur1 = mn + mur
    mnur2 = mn - mur

    G1_bin = tool_box.set_bin(mg1r, mg_bin_num, 100)
    G1_hist_bin = gglensing_tool.set_bin(mg1r, 4000, 1.001, "log")
    NU1_hist_bin = gglensing_tool.set_bin(mnur1, 4000, 1.001, "log")

    img = Image_Plot(xpad=0.25, ypad=0.25)
    img.subplots(1, 2)

    result = gglensing_tool.find_shear_grid(mg1r, mnur1, G1_bin, G1_hist_bin, NU1_hist_bin,
                                                             chisq_gap=50, dg=10, fit_num=20, ax=img.axs[0][0])[:4]
    gh, gh_sig, coeff, asym = result
    # chisqs_min = coeff[0] - coeff[1] ** 2 / 4 / coeff[2]
    # if ij[1] == 0:
    #     total_result[ij[0], 0] = asym
    #     total_result[ij[0]+dilute_case, 0] = chisqs_min

    print("%.3f(%.3f). %d source + %d contaminations" % (gh, gh_sig, num_s, num_non))

    # correction
    h5f = h5py.File(data_path3, "r")
    data = h5f["/data"][()]
    mg1, mg2, mn, mu, mv = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]
    h5f.close()
    h5f = h5py.File(para_path3, "r")
    ra = h5f["/ra"][()]
    dec = h5f["/dec"][()]
    z = h5f["/z"][()] + foreground_z_err
    h5f.close()

    com_dist = cosmos.comoving_distance(z).value * H_0 / 100

    crit_coeff = 1662895.2081868195 * com_dist / com_dist_len / (com_dist - com_dist_len) / (1 + len_z)

    src_pos = SkyCoord(ra=ra * units.arcsec, dec=dec * units.arcsec, frame="fk5")
    position_theta = len_pos.position_angle(src_pos).rad


    cos_2theta = numpy.cos(2 * position_theta)
    sin_2theta = numpy.sin(2 * position_theta)
    cos_4theta = numpy.cos(4 * position_theta)
    sin_4theta = numpy.sin(4 * position_theta)

    mg1r_corr = (mg1 * cos_2theta - mg2 * sin_2theta) * crit_coeff
    mg2r_corr = (mg1 * sin_2theta + mg2 * cos_2theta) * crit_coeff
    mur_corr = mu * cos_4theta - mv * sin_4theta
    mnur1_corr = mn + mur_corr
    mnur2_corr = mn - mur_corr


    result_c = gglensing_tool.find_shear_grid_corr_new(mg1r,
                                                 mnur1,
                                                 mg1r_corr[:num_corr],
                                                 mnur1_corr[:num_corr],
                                                 G1_bin, G1_hist_bin, NU1_hist_bin, chisq_gap=50, dg=10, fit_num=20,
                                                 ax=img.axs[0][1])[:4]
    ghc, ghc_sig, coeffc, asymc = result_c
    chisqs_min = coeffc[0] - coeffc[1] ** 2 / 4 / coeffc[2]
    total_result[ij[0], ij[1]] = asymc
    total_result[ij[0]+dilute_case, ij[1]] = chisqs_min

    print("%.3f(%.3f). %d source + %d contaminations + %d corrections." % (gh, gh_sig, num_s, num_non, num_corr))
    img.save_img("./%d/%d_source_%d_dilution_%d_corr.png"%(mg_bin_num, num_s, num_non, num_corr))
    # img.show_img()
    img.close_img()

comm.Barrier()

if rank == 0:
    numpy.savez("./%d/cache.npz"%mg_bin_num, total_result)

    img = Image_Plot(xpad=0.25,ypad=0.25)
    img.subplots(2,dilute_case)

    for i in range(dilute_case):
        img.axs[0][i].plot(corr_num, total_result[i])
        img.axs[1][i].plot(corr_num, total_result[i+dilute_case])

        xs = img.axs[0][i].set_xlim()
        img.axs[0][i].plot([xs[0],xs[1]], [total_result[i,0], total_result[i,0]],label="Asym before corr")

        xs = img.axs[1][i].set_xlim()
        img.axs[1][i].plot([xs[0],xs[1]], [total_result[i+dilute_case,0], total_result[i+dilute_case,0]],
                           label="$\chi^2$ before corr")


        ys = img.axs[0][i].set_ylim()
        dilute_ratio = dilute_ratio_list[i]
        num_non = int(num_s * dilute_ratio)
        img.axs[0][i].plot([num_non,num_non], [ys[0], ys[1]],ls="--",c="gray", label="true dilution")

        img.set_label(0,i,0,"Asymmetry")
        img.set_label(1,i,0,"$\chi^2$")

        img.set_label(0,i,1,"Correction")
        img.set_label(1,i,1,"Correction")

        img.axis_sci_ticklabel(0,i,0)
        img.axis_sci_ticklabel(0,i,1)

        img.axis_sci_ticklabel(1,i,1)

        img.axs[0][i].legend()
        img.axs[1][i].legend()


    img.save_img("./%d/asym_%d.pdf"%(mg_bin_num,mg_bin_num))
    img.close_img()
comm.Barrier()