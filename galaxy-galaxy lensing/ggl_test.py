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
from astropy.coordinates import SkyCoord
from astropy import units
import time
from mpi4py import MPI


def mix_data(rank, sheared_path, non_sheared_path, src_z_threshold, zerr_sig, rng, max_dilution=0.2):
    h5f_s = h5py.File(sheared_path, "r")
    src_z_s = h5f_s["/%d/z" % rank][()]

    if zerr_sig > 0.0001:
        src_z_err = rng.normal(0, (1 + src_z_s) * zerr_sig)
    else:
        src_z_err = 0

    src_z_ny = src_z_s + src_z_err

    idx = src_z_ny >= src_z_threshold

    src_num = idx.sum()
    dilution_num_max = int(src_num / (1 - max_dilution) * max_dilution)
    total_num = src_num + dilution_num_max

    src_z = numpy.zeros((total_num,), dtype=numpy.float32)
    src_ra = numpy.zeros((total_num,), dtype=numpy.float32)
    src_dec = numpy.zeros((total_num,), dtype=numpy.float32)
    src_radius = numpy.zeros((total_num, ), dtype=numpy.float32)
    src_radian = numpy.zeros((total_num, ), dtype=numpy.float32)
    mgt = numpy.zeros((total_num,), dtype=numpy.float32)
    mnu1 = numpy.zeros((total_num,), dtype=numpy.float32)


    src_z[:src_num] = src_z_ny[idx]
    src_ra[:src_num] = h5f_s["/%d/ra" % rank][()][idx]
    src_dec[:src_num] = h5f_s["/%d/dec" % rank][()][idx]
    src_radius[:src_num] = h5f_s["/%d/sep_radius" % rank][()][idx]
    src_radian[:src_num] = h5f_s["/%d/sep_radian" % rank][()][idx]
    mgt[:src_num] = h5f_s["/%d/mgt" % rank][()][idx]
    mnu1[:src_num] = h5f_s["/%d/mnu1" % rank][()][idx]
    h5f_s.close()

    # non-sheared data
    h5f_n = h5py.File(non_sheared_path, "r")

    src_ra[src_num:] = h5f_n["/%d/ra" % rank][()][:dilution_num_max]
    src_dec[src_num:] = h5f_n["/%d/dec" % rank][()][:dilution_num_max]
    src_radius[src_num:] = h5f_n["/%d/sep_radius" % rank][()][:dilution_num_max]
    src_radian[src_num:] = h5f_n["/%d/sep_radian" % rank][()][:dilution_num_max]

    z_n = numpy.abs(rng.normal(0, 0.1, 2*dilution_num_max)) + src_z_threshold
    src_z[src_num:] = z_n[:dilution_num_max]

    mgt[src_num:] = h5f_n["/%d/mgt" % rank][()][:dilution_num_max]
    mnu1[src_num:] = h5f_n["/%d/mnu1" % rank][()][:dilution_num_max]
    h5f_n.close()

    return src_z, src_ra, src_dec, src_radius,src_radian, mgt, mnu1, src_num




#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


data_path = argv[1]

cmd = int(argv[2])
chi_gap = 50

# zerr tag
zerr_sig = float(argv[3])

# "noise_free" or "noisy_cpp"
data_type = argv[4]

#"/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_1"
result_path = data_path + "/result/dilution_test/dilution_zerr_%.2f/%d"%(zerr_sig, cmd)

# bin number for PDF_SYM
pdf_bin_num = [2, 10]

dilution_ratio = numpy.array([0, 0.05, 0.1, 0.15, 0.2])

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

nfw_model = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)

# redshift threshold
src_z_threshold = len_z + 0.1



guess_num_p = 50
signal_guess = numpy.zeros((guess_num_p+guess_num_p,))
signal_guess[:guess_num_p] = -hk_tool_box.set_bin_log(0.001, 500, guess_num_p)
signal_guess[guess_num_p:] = hk_tool_box.set_bin_log(0.001, 500, guess_num_p)
signal_guess = numpy.sort(signal_guess)


rng = numpy.random.RandomState(213 + rank*124212)
# read the source data
sheared_path = data_path + "/data/segment_sheared_%s.hdf5"%data_type
non_sheared_path = data_path + "/data/segment_non_sheared_%s.hdf5"%data_type
# ra, dec Degree
src_z, src_ra, src_dec, src_radius, src_radian, mgt, mnu1, src_num = mix_data(rank,sheared_path, non_sheared_path, src_z_threshold, zerr_sig, rng)

total_num = (src_num/(1 - dilution_ratio)).astype(dtype=numpy.intc)

com_dist_src = cosmos.comoving_distance(src_z).value * h  # Mpc/h

mean_radius = src_radius.mean()
mean_radian = src_radian.mean()/numpy.pi*180*3600

# crit_sd_num = 1662895.2081868195*com_dist_src
# crit_sd_denorm = com_dist_len*(com_dist_src-com_dist_len)*(1+len_z)
# crit_sd = crit_sd_num/crit_sd_denorm
# crit_sd = 1662895.2081868195*com_dist_src/com_dist_len/(com_dist_src-com_dist_len)/(1+len_z)

sd_coeff = 1662895.2081868195
if cmd == 0:
    # coeff_1 = crit_sd_num/crit_sd_denorm
    coeff_1 = sd_coeff*com_dist_src/com_dist_len/(com_dist_src-com_dist_len)/(1+len_z)
    coeff_2 = 1
elif cmd == 1:
    coeff_1 = 1
    # coeff_2 = crit_sd_denorm/crit_sd_num
    coeff_2 = com_dist_len*(com_dist_src-com_dist_len)*(1+len_z)/sd_coeff/com_dist_src
else:
    # coeff_1 = crit_sd_num
    # coeff_2 = crit_sd_denorm
    coeff_1 = sd_coeff*com_dist_src
    coeff_2 = com_dist_len*(com_dist_src-com_dist_len)*(1+len_z)


mgt = mgt*coeff_1
mnu1 = mnu1*coeff_2


#  Model
Ds_true = hk_gglensing_tool.get_delta_sigma(nfw_model, com_dist_len, len_z, com_dist_src[0], src_z[0], numpy.array([mean_radian]))


t1 = time.time()


result_sp = [len(pdf_bin_num), 5*len(dilution_ratio)]
result_sub = numpy.zeros((result_sp[0], result_sp[1]))


chisq_img = Image_Plot(ypad=0.25,xpad=0.2)
chisq_img.subplots(len(pdf_bin_num), len(dilution_ratio))

for i in range(len(dilution_ratio)):

    st, ed = int(i * 5), int((i + 1) * 5)
    temp = numpy.random.choice(mgt[:total_num[i]], 100000, False)

    print("Totla num: %d. src num: %d(%.2f). "
          "Dilution: %d(%.2f)."% (total_num[i], src_num, src_num/total_num[i], total_num[i]-src_num, (total_num[i]-src_num)/total_num[i]))
    for j in range(len(pdf_bin_num)):
        # result_t = hk_FQlib.find_shear_cpp(mgt_nf[idx], mnu1_nf[idx], bin_num=pdf_bin_num[j], left=-200, right=200,
        #                                    chi_gap=chi_gap,max_iters=60,fig_ax=chisq_img.axs[0][j])[:2]

        pdf_bin = hk_FQlib.set_bin_(temp, pdf_bin_num[j], scale=1000000)

        result_t = hk_FQlib.find_shear_cpp_guess(mgt[:total_num[i]], mnu1[:total_num[i]],
                                                 pdf_bin, signal_guess, chi_gap=chi_gap, fig_ax=chisq_img.axs[j][i])

        result_sub[j, st:ed] = total_num[i],mean_radius, Ds_true[0], result_t[0], result_t[1]
        chisq_img.axs[j][i].set_title("dilution_ratio:%.2f. %d PDF_bin. True signal: %.4f"%(dilution_ratio[i], pdf_bin_num[j], Ds_true[0]))

chisq_img.save_img(result_path + "/imgs/%d_chisq_%s.pdf"%(rank, data_type))
chisq_img.close_img()

t2 = time.time()
comm.Barrier()
print(rank, t2-t1)

if rank > 0:
    # !!!! remember the data type, MPI.DOUBLE, MPI.FLOAT, ...
    # or it will raise an error, Keyerror
    comm.Send([result_sub, MPI.DOUBLE], dest=0, tag=rank)
else:

    result_collect = numpy.zeros((numprocs, result_sp[0], result_sp[1]))
    result_collect[0, :, :] = result_sub

    # receive the data from other CPUs
    # !!!! the start points is 1 in range() not 0
    for procs in range(1, numprocs):
        # prepare a buffer for the data, the shape must be the same
        # with that of what the other CPUs send, you have collected them in 'data_sps'
        recvs = numpy.empty((result_sp[0], result_sp[1]), dtype=numpy.double)
        # receive it using the buffer,
        comm.Recv(recvs, source=procs, tag=procs)

        result_collect[procs,:,:] = recvs

    numpy.savez(result_path + "/result_%s.npz"%data_type, result_collect)

    # img = Image_Plot(xpad=0.25, ypad=0.24)
    # img.subplots(len(pdf_bin_num), len(dilution_ratio))
    #
    # for j in range(len(pdf_bin_num)):
    #     tag = int(4 * j)
    #     img.axs[j][0].errorbar(inform[0], Ds_nf[tag], Ds_nf[tag + 1], capsize=3, marker="s", fmt=" ",
    #                            label="Noise free")
    #     img.axs[j][0].errorbar(inform[0], Ds_ny[tag], Ds_ny[tag + 1], capsize=3, marker="o", fmt=" ",
    #                            label="Noisy")
    #
    #     img.axs[j][1].errorbar(inform[0], Ds_nf[tag + 2], Ds_nf[tag + 3], capsize=3, marker="s", fmt=" ",
    #                            label="Noise free")
    #     img.axs[j][1].errorbar(inform[0], Ds_ny[tag + 2], Ds_ny[tag + 3], capsize=3, marker="o", fmt=" ",
    #                            label="Noisy")
    #
    #     img.axs[j][0].plot(inform[0], Ds_true, label="Model")
    #
    #     img.axs[j][2].plot(inform[0], (Ds_nf[tag] - Ds_true), marker="s", label="Noise free")
    #     img.axs[j][2].plot(inform[0], (Ds_ny[tag] - Ds_true), marker="o",label="Noisy")
    #     # img.axs[0][2].set_yscale("symlog")
    #
    #     img.axs[j][3].plot(inform[0], (Ds_nf[tag] - Ds_true) / Ds_nf[tag+1],marker="s", label="Noise free")
    #     img.axs[j][3].plot(inform[0], (Ds_ny[tag] - Ds_true) / Ds_ny[tag+1], marker="o",label="Noisy")
    #
    #
    #     # img.axs[j][2].plot(inform[0], Ds_nf[tag], marker="s", label="Noise free")
    #     # img.axs[j][2].plot(inform[0], Ds_ny[tag], marker="o",label="Noisy")
    #     # # img.axs[0][2].set_yscale("symlog")
    #     #
    #     # img.axs[j][3].plot(inform[0], Ds_nf[tag] / Ds_nf[tag+1],marker="s", label="Noise free")
    #     # img.axs[j][3].plot(inform[0], Ds_ny[tag] / Ds_ny[tag+1], marker="o",label="Noisy")
    #
    #
    #     for i in range(4):
    #         img.axs[j][i].set_title("PDF bin_num %d"%pdf_bin_num[j])
    #         img.set_label(j, i, 1, "Radius Mpc/h")
    #         img.axs[j][i].set_xscale("log")
    #         img.axs[j][i].legend()
    #     img.set_label(j, 0, 0, "$\Delta\Sigma$")
    #     img.set_label(j, 1, 0, "$\Delta\Sigma_x$")
    #
    #     img.set_label(j, 2, 0, "$\Delta\Sigma - \Delta\Sigma_{model}$")
    #     img.set_label(j, 3, 0, "$(\Delta\Sigma - \Delta\Sigma_{model})/Error bar$")
    #
    #     img.axs[j][0].set_yscale("log")
    # img.save_img(result_path + "/signal_comparison_%d_%d.pdf"%(numprocs,tt))
    # img.show_img()

    # img = Image_Plot(xpad=0.25)
    # img.subplots(1, 2)
    # for i in range(1, 3):
    #     tag = int(i * 4)
    #     img.axs[0][0].plot(inform[0], Ds_nf[tag + 1] / Ds_nf[1],
    #                        label="Noise free: Error bar %d bins/ %d bins" % (pdf_bin_num[i], pdf_bin_num[0]))
    #     img.axs[0][1].plot(inform[0], Ds_ny[tag + 1] / Ds_ny[1],
    #                        label="Noisy: Error bar %d bins/ %d bins" % (pdf_bin_num[i], pdf_bin_num[0]))
    # for i in range(2):
    #     img.axs[0][i].legend()
    #     img.axs[0][i].set_xscale("log")
    #     img.set_label(0, i, 0, "Error bar ratio")
    #     img.set_label(0, i, 1, "Radius Mpc/h")
    # img.save_img(result_path + "/err_comparison_%d.pdf"%numprocs)
    # img.show_img()

    # img = Image_Plot()
    # img.subplots(1,1)
    # img.axs[0][0].scatter(src_z_true[:1000],src_z[:1000])
    # x1,x2 = min(src_z_true[:1000].min(), src_z[:1000].min()),max(src_z_true[:1000].max(), src_z[:1000].max())
    #
    # img.axs[0][0].plot([x1,x2],[x1,x2],ls="dashed",c="k")
    # img.set_label(0,0,0,"Z")
    # img.set_label(0,0,1,"Z_true")
    # img.save_img(result_path + "/z.pdf")

comm.Barrier()

