from sys import path, argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import gglensing_tool
import numpy
import h5py
import galsim
import FQlib
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
import time


def get_log(x):
    idx1 = x < 0
    idx2 = x > 0
    logx = numpy.zeros_like(x)
    logx[idx1] = -numpy.log(-x[idx1])
    logx[idx2] = numpy.log(x[idx2])
    return logx


def get_diff_1d(hist1d):
    hist_num = hist1d.shape[0]
    hist_num2 = int(hist_num / 2)

    hist1d_left = numpy.flip(hist1d[:hist_num2], axis=0)
    hist1d_right = hist1d[hist_num2:]

    hist1d_diff = hist1d_right - hist1d_left
    hist1d_sum = hist1d_right + hist1d_left

    idx = hist1d_sum > 0
    chisq = numpy.zeros_like(hist1d_left)
    chisq[idx] = hist1d_diff[idx] ** 2 / hist1d_sum[idx] / 2

    return chisq, hist1d_diff, hist1d_sum, hist1d_left, hist1d_right


def get_diff(hist2d):
    hist_num = hist2d.shape[1]
    hist_num2 = int(hist_num / 2)

    hist2d_left = numpy.flip(hist2d[:, :hist_num2], axis=1)
    hist2d_right = hist2d[:, hist_num2:]

    hist2d_diff = hist2d_right - hist2d_left
    hist2d_sum = hist2d_right + hist2d_left

    idx = hist2d_sum > 0
    chisq = numpy.zeros_like(hist2d_left)
    chisq[idx] = hist2d_diff[idx] ** 2 / hist2d_sum[idx] / 2

    return chisq, hist2d_left, hist2d_right, hist2d_diff, hist2d_sum


def prepare_data(data_path, para_path, cosmos, len_z, H_0, len_pos, dz=False, zerr=0, scale=1.):
    h5f = h5py.File(data_path, "r")
    data = h5f["/data"][()]*scale
    h5f.close()

    h5f = h5py.File(para_path, "r")
    ra = h5f["/ra"][()]
    dec = h5f["/dec"][()]

    if dz:
        z = zerr
    else:
        z = h5f["/z"][()]

    h5f.close()

    com_dist = cosmos.comoving_distance(z).value * H_0 / 100

    com_dist_len = cosmos.comoving_distance(len_z).value * H_0 / 100

    crit_coeff = 1662895.2081868195 * com_dist / com_dist_len / (com_dist - com_dist_len) / (1 + len_z)

    src_pos = SkyCoord(ra=ra * units.arcsec, dec=dec * units.arcsec, frame="fk5")
    position_theta = len_pos.position_angle(src_pos).rad

    mg1r, mg2r, mnur1, mnur2 = FQlib.rotate(data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4], position_theta)
    mg1r *= crit_coeff
    print(mnur1.min(), mnur1.max())

    # the true signal
    separation_theta = len_pos.separation(src_pos).arcsec.mean()
    zm = z.mean()
    com_dist_m = cosmos.comoving_distance(zm).value * H_0 / 100
    delta_sigma = gglensing_tool.get_delta_sigma(nfw, com_dist_len, len_z, com_dist_m, zm, numpy.array([separation_theta]))

    return mg1r, mnur1, crit_coeff, delta_sigma[0]


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
print("Lens plane at z = %.2f, %.5f Mpc/h"%(len_z, com_dist_len))
len_pos = SkyCoord(ra=0*units.arcsec, dec=0*units.arcsec,frame="fk5")

# lens profile
nfw = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)


# data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z"
data_path = "/coma/hklee/Galaxy_Galaxy_lensing_test/background/continue_source_z"

hist_bin_tag = 0
shear_tag = 0

rng = numpy.random.RandomState(12411)
foreground_z_err = numpy.abs(rng.normal(0, 0.2, 20000000)) + 0.31
# foreground_z_err = numpy.abs(rng.normal(0, 0.5, 20000000)) + 0.31


num_total = 20000000
num_src = 16000000
num_non = num_total - num_src

dilution_ratio = num_non*1.0/num_total

mg_pdf_bin_num = 50
mg_pdf_bin_num2 = int(mg_pdf_bin_num/2)

data_path1 = data_path + "/data/sheared_data_%d_noise_free.hdf5"%shear_tag
data_path2 = data_path + "/data/non_sheared_data_%d_noise_free.hdf5"%shear_tag
data_path3 = data_path + "/data/non_sheared_data_%d_noise_free.hdf5"%(1+shear_tag)

para_path1 = data_path + "/paras/sheared_para_%d.hdf5"%shear_tag
para_path2 = data_path + "/paras/non_sheared_para_%d.hdf5"%shear_tag
para_path3 = data_path + "/paras/non_sheared_para_%d.hdf5"%(1+shear_tag)

mg1r_src, mnur1_src, crit_coeff_src, delta_sigma_src = prepare_data(data_path1, para_path1, cosmos, len_z, H_0,
                                                                    len_pos)

mg1r_non, mnur1_non, crit_coeff_non = prepare_data(data_path2, para_path2, cosmos, len_z, H_0,
                                                    len_pos, True, foreground_z_err, 1)[:3]

mg1r_corr_1, mnur1_corr_1, crit_coeff_corr_1 = prepare_data(data_path3, para_path3, cosmos, len_z, H_0,
                                                            len_pos, True,  foreground_z_err, 0.5)[:3]

mg1r_corr_2, mnur1_corr_2, crit_coeff_corr_2 = prepare_data(data_path3, para_path3, cosmos, len_z, H_0,
                                                            len_pos, True, foreground_z_err, 1)[:3]

print("The true \Delta\Simga=%.4f"%delta_sigma_src)

mg_pdf_bin_1d = FQlib.set_bin_(mg1r_src, mg_pdf_bin_num, 100)


mg1r_total = numpy.zeros((num_total, ))
mnur1_total = numpy.zeros((num_total, ))

mg1r_total[:num_src] = mg1r_src[:num_src]
mg1r_total[num_src:] = mg1r_non[:num_non]

mnur1_total[:num_src] = mnur1_src[:num_src]
mnur1_total[num_src:] = mnur1_non[:num_non]


# set bins for 2d histogram
hist2d_xbin_num = 200
hist2d_xbin_num2 = int(hist2d_xbin_num/2)

if mnur1_total.min() < 0:
    hist2d_ybin_num = 200
else:
    hist2d_ybin_num = 100

mg_bin = gglensing_tool.set_bin(mg1r_src, hist2d_xbin_num, bound_scale=1.1, method="log")
mnur_bin = gglensing_tool.set_bin(mnur1_src, hist2d_ybin_num, bound_scale=1.1, method="log")

x = numpy.zeros((hist2d_ybin_num,hist2d_xbin_num))
y = numpy.zeros((hist2d_ybin_num,hist2d_xbin_num))

for i in range(hist2d_xbin_num):
    x[:,i] = (mg_bin[i] + mg_bin[i+1])/2
for i in range(hist2d_ybin_num):
    y[i] = (mnur_bin[i] + mnur_bin[i+1])/2

xh, yh = x[:,hist2d_xbin_num2:], y[:,hist2d_xbin_num2:]

xhlog, yhlog = get_log(xh), get_log(yh)
xlog, ylog = get_log(x), get_log(y)



# # plot the distribution of \Sigma_{crit} and G_t
# img = Image_Plot(xpad=0.25, ypad=0.25)
# img.subplots(2,2)
#
# mg1r_src_hist = numpy.histogram(mg1r_src, mg_bin)[0]
# mg1r_non_hist = numpy.histogram(mg1r_non, mg_bin)[0]
# idx1 = mg1r_src_hist > 0
# idx2 = mg1r_non_hist > 0
# idx = idx1 | idx2
# img.axs[0][0].scatter(x[0][idx1], mg1r_src_hist[idx1], s=5, label="background $\Sigma_{crit}G_t$")
# img.axs[0][0].scatter(x[0][idx2], mg1r_non_hist[idx2], s=5,label="contamination $\Sigma_{crit}G_t$")
#
# img.axs[1][0].scatter(x[0][idx], mg1r_src_hist[idx] - mg1r_non_hist[idx],s=5,label="diff $\Sigma_{crit}G_t$")
#
#
# hist_crit_src, crit_bins = numpy.histogram(crit_coeff_src, 2000)[:2]
# hist_crit_non = numpy.histogram(crit_coeff_non, crit_bins)[0]
#
# crit_bins_cent = (crit_bins[1:] + crit_bins[:-1])/2
#
# img.axs[0][1].scatter(crit_bins_cent, hist_crit_src, s=5,label="background $\Sigma_{crit}$")
# img.axs[0][1].scatter(crit_bins_cent, hist_crit_non,s=5, label="contamination $\Sigma_{crit}$")
# img.axs[1][1].scatter(crit_bins_cent, hist_crit_src - hist_crit_non, s=5,label="diff $\Sigma_{crit}$")
#
# for i in range(2):
#     # img.axs[i][0].set_xscale("symlog")
#     img.axs[0][i].set_yscale("log")
#     # img.axs[1][i].set_yscale("symlog")
#     for j in range(2):
#         img.axs[i][j].legend()
#
# img.save_img("./crit_hist.png")
#
#
img = Image_Plot(xpad=0.25, ypad=0.25)
img.subplots(1, 5)

t1 = time.time()

src_result = FQlib.find_shear_cpp(mg1r_src, mnur1_src, mg_pdf_bin_num, left=-100, right=200,
                                  fit_num=30, chi_gap=50,fig_ax=img.axs[0][0])
src_Ds, src_Ds_err, src_coeffs, src_chisqs_min, src_bins = src_result

t2 = time.time()

print("%d foreground. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_src,src_Ds, src_Ds_err, src_chisqs_min, t2-t1))

non_result = FQlib.find_shear_cpp(mg1r_non, mnur1_non, mg_pdf_bin_num, left=-100, right=200,
                                  fit_num=30, chi_gap=50,fig_ax=img.axs[0][1])
non_Ds, non_Ds_err,non_coeffs, non_chisqs_min, non_bins = non_result

t3 = time.time()
print("%d foreground. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_src, non_Ds, non_Ds_err, non_chisqs_min, t3-t2))

total_result = FQlib.find_shear_cpp(mg1r_total, mnur1_total, mg_pdf_bin_num, left=-100, right=200,
                                  fit_num=30, chi_gap=50,fig_ax=img.axs[0][2])
total_Ds, total_Ds_err,total_coeffs, total_chisqs_min, total_bins = total_result

t4 = time.time()
print("%d foreground + %d contamination. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_src, num_non, total_Ds,
                                                                                 total_Ds_err, total_chisqs_min, t4-t3))


# correction using similar G and N

total_result_corr_1 = FQlib.find_shear_cpp_corr(mg1r_total, mnur1_total, mg1r_corr_1[:num_non], mnur1_corr_1[:num_non],
                                                mg_pdf_bin_num, left=-100, right=200, fit_num=30, chi_gap=50,fig_ax=img.axs[0][3])
total_Ds_corr_1, total_Ds_err_corr_1, total_coeffs_corr_1, total_chisqs_min_corr_1, total_bins_corr_1 = total_result_corr_1

t5 = time.time()
print("%d foreground + %d contamination. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_src, num_non, total_Ds_corr_1,
                                                                                 total_Ds_err_corr_1, total_chisqs_min_corr_1, t5-t4))


# correction using different G and N

total_result_corr_2 = FQlib.find_shear_cpp_corr(mg1r_total, mnur1_total, mg1r_corr_2[:num_non], mnur1_corr_2[:num_non],
                                                mg_pdf_bin_num, left=-100, right=200, fit_num=30, chi_gap=50,fig_ax=img.axs[0][4])
total_Ds_corr_2, total_Ds_err_corr_2,total_coeffs_corr_2, total_chisqs_min_corr_2, total_bins_corr_2 = total_result_corr_2

t6 = time.time()
print("%d foreground + %d contamination. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_src, num_non, total_Ds_corr_2,
                                                                                 total_Ds_err_corr_2, total_chisqs_min_corr_2, t6-t5))


for i in range(5):
    img.set_label(0,i,0,"$\chi^2$")
    img.set_label(0,i,1,"$\Delta\Sigma$")
img.save_img("./pdf_sym.png")
img.show_img()
img.close_img()



# # plot the count diff and chi^2 before and after PDF_SYM
#
# mg_pdf_bin_label = [-i for i in range(1,mg_pdf_bin_num2+1)]
# for i in range(1,mg_pdf_bin_num2+1):
#     mg_pdf_bin_label.append(i)
# mg_pdf_bin_label = numpy.sort(numpy.array(mg_pdf_bin_label))
#
# img = Image_Plot(xpad=0.25, ypad=0.25)
# img.subplots(2,3)
# total_Ds = 59.491
# for i in range(2):
#
#     corr_mg1r = mg1r_total - i*total_Ds*mnur1_total
#
#     hist_total = numpy.histogram(corr_mg1r, mg_pdf_bin_1d)[0]
#     chisq_total, hist_diff_total = get_diff_1d(hist_total)[:2]
#
#     hist_src = numpy.histogram(corr_mg1r[:num_src], mg_pdf_bin_1d)[0]
#     chisq_src, hist_diff_src = get_diff_1d(hist_src)[:2]
#
#     hist_non = numpy.histogram(corr_mg1r[num_src:], mg_pdf_bin_1d)[0]
#     chisq_non, hist_diff_non = get_diff_1d(hist_non)[:2]
#
#     img.axs[i][0].scatter(mg_pdf_bin_label, hist_total, label="Total. hist")
#     img.axs[i][0].scatter(mg_pdf_bin_label, hist_src, label="back. hist")
#     img.axs[i][0].scatter(mg_pdf_bin_label, hist_non, label="fore. hist")
#
#
#     img.axs[i][1].scatter(numpy.arange(1, mg_pdf_bin_num2+1), hist_diff_total, label="Total. $n_{right} - n_{left}$")
#     img.axs[i][1].scatter(numpy.arange(1, mg_pdf_bin_num2+1), hist_diff_src, label="back. $n_{right} - n_{left}$")
#     img.axs[i][1].scatter(numpy.arange(1, mg_pdf_bin_num2+1), hist_diff_non, label="fore. $n_{right} - n_{left}$")
#
#     img.axs[i][2].scatter(numpy.arange(1, mg_pdf_bin_num2+1), chisq_total, label="Total. $\chi^2$")
#     img.axs[i][2].scatter(numpy.arange(1, mg_pdf_bin_num2+1), chisq_src, label="back. $\chi^2$")
#     img.axs[i][2].scatter(numpy.arange(1, mg_pdf_bin_num2+1), chisq_non, label="fore. $\chi^2$")
#
#     for j in range(3):
#         if i == 0:
#             img.axs[i][j].set_title("before PDF_SYM")
#         else:
#             img.axs[i][j].set_title("after PDF_SYM")
#         img.set_label(i,j, 1, "bin label")
#         img.axs[i][j].legend()
#         img.axis_sci_ticklabel(i, j, 0)
#     img.set_label(i, 1, 0, "$n_{right} - n_{left}$")
#     img.set_label(i, 2, 0, "$\chi^2$")
#
# img.save_img("./components_chisq.png")
#
# exit()

cmd = argv[1]
if cmd == "src":
    check_mg = mg1r_src
    check_mnu = mnur1_src
    # sigma = 72.816
    sigma = delta_sigma_src

elif cmd == "non":
    check_mg = mg1r_non
    check_mnu = mnur1_non
    sigma = -0.198
else:
    check_mg = mg1r_total
    check_mnu = mnur1_total
    sigma = 56.48

print(cmd, sigma, delta_sigma_src)

# hist_bin_label = [-i for i in range(1,hist2d_bin_num2+1)]
# for i in range(1,hist2d_bin_num2+1):
#     hist_bin_label.append(i)
# hist_bin_label = numpy.sort(numpy.array(hist_bin_label))
#
# img = Image_Plot(xpad=0.25, ypad=0.25)
# img.subplots(1, 3)
# for i in range(2):
#
#     hist = numpy.histogram(check_mg - i*sigma*check_mnu, mg_bin)[0]
#     n1, n2 = numpy.flip(hist[:hist2d_bin_num2],axis=0), hist[hist2d_bin_num2:]
#
#     hist_diff = n2 - n1
#     hist_sum = n2 + n1
#
#     idx = hist_sum > 0
#     chisq = hist_diff[idx]**2/hist_sum[idx]
#
#     idx1 = hist > 0
#     img.axs[0][0].scatter(hist_bin_label[idx1], hist[idx1], label="$\chi^2=%.3f$"%numpy.sum(chisq))
#
#     img.axs[0][1].scatter(hist_bin_label[hist2d_bin_num2:][idx], hist_diff[idx],label="$\chi^2=%.3f$"%numpy.sum(chisq))
#     img.axs[0][2].scatter(hist_bin_label[hist2d_bin_num2:][idx], chisq,label="$\chi^2=%.3f$"%numpy.sum(chisq))
#
# for i in range(3):
#
#     if i == 2:
#         ys = img.axs[0][i].set_ylim()
#         img.axs[0][i].set_ylim((0.1,ys[1]))
#         img.axs[0][i].set_yscale("log")
#     img.axs[0][i].legend()
#
# img.save_img("./chisq_src.png")


img = Image_Plot(xpad=0.25, ypad=0.25)
img.subplots(2, 3)

titles = [["Hist before PDF_SYM", "$\Delta N$ before PDF_SYM", "$\chi^2$ before PDF_SYM"],
            ["Hist after PDF_SYM", "$\Delta N$ after PDF_SYM", "$\chi^2$ after PDF_SYM"]]


for i in range(2):

    hist2d = numpy.histogram2d(check_mnu, check_mg - i*sigma*check_mnu, [mnur_bin, mg_bin])[0]
    chisq, hist2d_left, hist2d_right, hist2d_diff, hist2d_sum = get_diff(hist2d)

    numpy.savez("./cache_%d_%s.npz"%(i, cmd), x,y, xh, yh, hist2d, chisq, hist2d_left, hist2d_right, hist2d_diff, hist2d_sum)

    idx = hist2d > 0
    img.scatter_pts(i, 0, xlog[idx], xlog[idx], hist2d[idx], color_map="jet")

    idx = hist2d_sum > 0
    img.scatter_pts(i, 1, xhlog[idx], yhlog[idx], hist2d_diff[idx], color_map="jet")

    img.scatter_pts(i, 2, xhlog[idx], yhlog[idx], chisq[idx], color_map="jet")

    for j in range(3):
        img.set_label(i,j,0,"N+U bin [log]")
        img.set_label(i,j,1,"$G_t$ bin [log]")
        img.axs[i][j].set_title(titles[i][j])
img.save_img("./pdf_sym_hist2d_%s.png"%cmd)


img = Image_Plot(xpad=0.25, ypad=0.25)
img.subplots(2, 3)

for i in range(2):

    hist2d = numpy.histogram2d(check_mnu, check_mg - i*sigma*check_mnu, [mnur_bin, mg_bin])[0]
    chisq, hist2d_left, hist2d_right, hist2d_diff, hist2d_sum = get_diff(hist2d)

    idx = hist2d <= 0
    hist2d[idx] = numpy.nan
    fig = img.axs[i][0].imshow(numpy.flip(hist2d, axis=0), cmap="jet")
    img.figure.colorbar(fig,ax=img.axs[i][0])

    idx = hist2d_sum <= 0
    hist2d_diff[idx] = numpy.nan
    fig = img.axs[i][1].imshow(numpy.flip(hist2d_diff, axis=0), cmap="jet")
    img.figure.colorbar(fig,ax=img.axs[i][1])

    chisq[idx] = numpy.nan
    fig = img.axs[i][2].imshow(numpy.flip(chisq, axis=0), cmap="jet")
    img.figure.colorbar(fig,ax=img.axs[i][2])

    for j in range(3):
        img.set_label(i,j,0,"N+U bin")
        img.set_label(i,j,1,"$G_t$ bin")
        img.axs[i][j].set_title(titles[i][j])
img.save_img("./pdf_sym_arr_hist2d_%s.png"%cmd)