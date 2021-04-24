from sys import path, argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import gglensing_tool
import numpy
import h5py
import galsim
import FQlib
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
from mpi4py import MPI
import time





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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


itemsize = MPI.DOUBLE.Get_size()
element_num = 31
if rank == 0:
    # bytes for 10 double elements
    nbytes = 4*element_num*element_num*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)

buf1, itemsize = win1.Shared_query(0)

final_result = numpy.ndarray(buffer=buf1, dtype='d', shape=(int(4*element_num), element_num))



test_scale = numpy.linspace(0.2, 5, element_num)
scale_gridx = numpy.zeros((element_num, element_num))
scale_gridy = numpy.zeros((element_num, element_num))
for i in range(element_num):
    scale_gridy[:,i] = test_scale
    scale_gridx[i] = test_scale

scale_gridx = scale_gridx.flatten()
scale_gridy = scale_gridy.flatten()

test_pairs = [[i,j] for i in range(element_num) for j in range(element_num)]
pairs_sub = tool_box.alloc(test_pairs, numprocs)[rank]



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


data_type = argv[1]

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



data_path1 = data_path + "/data/sheared_data_%d_%s.hdf5"%(shear_tag, data_type)
data_path2 = data_path + "/data/non_sheared_data_%d_%s.hdf5"%(shear_tag, data_type)
data_path3 = data_path + "/data/non_sheared_data_%d_%s.hdf5"%(1+shear_tag,data_type)

para_path1 = data_path + "/paras/sheared_para_%d.hdf5"%(shear_tag)
para_path2 = data_path + "/paras/non_sheared_para_%d.hdf5"%(shear_tag)
para_path3 = data_path + "/paras/non_sheared_para_%d.hdf5"%(1+shear_tag)

mg1r_src, mnur1_src, crit_coeff_src, delta_sigma = prepare_data(data_path1, para_path1, cosmos, len_z, H_0, len_pos, False)

mg1r_non, mnur1_non, crit_coeff_non = prepare_data(data_path2, para_path2, cosmos, len_z, H_0,
                                                   len_pos, True, foreground_z_err, 1)[:3]

mg1r_corr, mnur1_corr, crit_coeff_corr = prepare_data(data_path3, para_path3, cosmos, len_z, H_0,
                                                            len_pos, True, foreground_z_err, 1)[:3]


mg1r_total = numpy.zeros((num_total,))
mnur1_total = numpy.zeros((num_total,))



mg_pdf_bin_num = 50
mg_pdf_bin_num2 = int(mg_pdf_bin_num/2)

mg_pdf_bin_1d = FQlib.set_bin_(mg1r_src, mg_pdf_bin_num, 100)

mg_pdf_bin_label = [-i for i in range(1,mg_pdf_bin_num2+1)]
for i in range(1,mg_pdf_bin_num2+1):
    mg_pdf_bin_label.append(i)
mg_pdf_bin_label = numpy.sort(numpy.array(mg_pdf_bin_label))




# set bins for 2d histogram
hist2d_bin_num = 200
hist2d_bin_num2 = int(hist2d_bin_num / 2)
mg_bin = gglensing_tool.set_bin(mg1r_src, hist2d_bin_num, bound_scale=1.1, method="log")
mnur_bin = gglensing_tool.set_bin(mnur1_src, hist2d_bin_num, bound_scale=1.1, method="log")

x = numpy.zeros((hist2d_bin_num, hist2d_bin_num))
y = numpy.zeros((hist2d_bin_num, hist2d_bin_num))

for i in range(hist2d_bin_num):
    x[:, i] = (mg_bin[i] + mg_bin[i + 1]) / 2
    y[i] = (mnur_bin[i] + mnur_bin[i + 1]) / 2

xh, yh = x[:, hist2d_bin_num2:], y[:, hist2d_bin_num2:]


for ij in pairs_sub:
    tag_i, tag_j = ij

    scale_i, scale_j = test_scale[tag_i], test_scale[tag_j]

    mg1r_total[:num_src] = mg1r_src[:num_src]
    mg1r_total[num_src:] = mg1r_non[:num_non]*scale_i

    mnur1_total[:num_src] = mnur1_src[:num_src]
    mnur1_total[num_src:] = mnur1_non[:num_non]*scale_i


    img = Image_Plot(xpad=0.25, ypad=0.25)
    img.subplots(1, 2)

    t3 = time.time()

    total_result = FQlib.find_shear_cpp(mg1r_total, mnur1_total, mg_pdf_bin_num, left=-100, right=200,
                                      fit_num=20, chi_gap=50,fig_ax=img.axs[0][0])
    total_Ds, total_Ds_err,total_coeffs, total_chisqs_min, total_bins = total_result

    t4 = time.time()
    print("%d foreground + %d contamination. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_src, num_non, total_Ds,
                                                                                     total_Ds_err, total_chisqs_min, t4-t3))


    # correction using similar G and N

    total_result_corr_1 = FQlib.find_shear_cpp_corr(mg1r_total, mnur1_total,
                                                    mg1r_corr[:num_non]*scale_j, mnur1_corr[:num_non]*scale_j,
                                                    mg_pdf_bin_num, left=-100, right=200, fit_num=20, chi_gap=50,
                                                    fig_ax=img.axs[0][1])
    total_Ds_corr_1, total_Ds_err_corr_1, total_coeffs_corr_1, total_chisqs_min_corr_1, total_bins_corr_1 = total_result_corr_1

    t5 = time.time()
    print("%d foreground + %d contamination. %.2f(%.2f). chisq_min: %.2f. %.2f sec"%(num_src, num_non, total_Ds_corr_1,
                                                                                     total_Ds_err_corr_1, total_chisqs_min_corr_1, t5-t4))

    final_result[tag_i,tag_j] = total_Ds
    final_result[tag_i + element_num, tag_j] = total_chisqs_min
    final_result[tag_i + int(2*element_num), tag_j] = total_Ds_corr_1
    final_result[tag_i + int(3*element_num), tag_j] = total_chisqs_min_corr_1

    for i in range(2):
        img.set_label(0,i,0,"$\chi^2$")
        img.set_label(0,i,1,"$\Delta\Sigma$")

        img.axs[0][i].set_title("%.4f-%.4f"%(scale_i, scale_j))

    img.save_img("./result/pdf/pdf_sym_%s_%d_%d.png"%(data_type, tag_i,tag_j))
    # img.show_img()
    img.close_img()

    # histogram before PDF_SYM
    hist1d_ori = numpy.histogram(mg1r_total, mg_pdf_bin_1d)[0]
    hist1d_chisq_ori, hist1d_diff_ori = get_diff_1d(hist1d_ori)[:2]

    hist2d_ori = numpy.histogram2d(mnur1_total, mg1r_total, [mnur_bin, mg_bin])[0]
    hist2d_chisq_ori, hist2d_diff_ori = get_diff(hist2d_ori)[:2]


    # histogram after PDF_SYM
    mg1r_sym = mg1r_total - mnur1_total*total_Ds

    hist1d_ori_sym = numpy.histogram(mg1r_sym, mg_pdf_bin_1d)[0]
    hist1d_chisq_ori_sym, hist1d_diff_ori_sym = get_diff_1d(hist1d_ori_sym)[:2]

    hist2d_ori_sym = numpy.histogram2d(mnur1_total, mg1r_sym, [mnur_bin, mg_bin])[0]
    hist2d_chisq_ori_sym, hist2d_diff_ori_sym = get_diff(hist2d_ori_sym)[:2]


    # corrected histogram after PDF_SYM
    mg1r_corr_sym1 = mg1r_total - mnur1_total*total_Ds_corr_1
    mg1r_corr_sym2 = mg1r_corr[:num_non]*scale_j - mnur1_corr[:num_non]*scale_j*total_Ds_corr_1

    hist1d_corr_sym1 = numpy.histogram(mg1r_corr_sym1, mg_pdf_bin_1d)[0]
    hist1d_corr_sym2 = numpy.histogram(mg1r_corr_sym2, mg_pdf_bin_1d)[0]
    hist1d_chisq_corr, hist1d_diff_corr = get_diff_1d(hist1d_corr_sym1 - hist1d_corr_sym2)[:2]

    hist2d_corr_sym1 = numpy.histogram2d(mnur1_total, mg1r_corr_sym1, [mnur_bin, mg_bin])[0]
    hist2d_corr_sym2 = numpy.histogram2d(mnur1_corr[:num_non]*scale_j, mg1r_corr_sym2, [mnur_bin, mg_bin])[0]
    hist2d_chisq_corr, hist2d_diff_corr = get_diff(hist2d_corr_sym1 - hist2d_corr_sym2)[:2]


    # # plot the count diff and chi^2 before and after PDF_SYM
    # img = Image_Plot(xpad=0.25, ypad=0.25)
    # img.subplots(2,3)
    #
    #
    # idx = hist2d_chisq_ori > 0
    # img.scatter_pts(0,0,xh[idx], yh[idx], hist2d_diff_ori[idx],color_map="jet")
    # img.scatter_pts(1,0,xh[idx], yh[idx], hist2d_chisq_ori[idx],color_map="jet")
    #
    #
    # idx = hist2d_chisq_ori_sym > 0
    # img.scatter_pts(0,1,xh[idx], yh[idx], hist2d_diff_ori_sym[idx],color_map="jet")
    # img.scatter_pts(1,1,xh[idx], yh[idx], hist2d_chisq_ori_sym[idx],color_map="jet")
    #
    #
    # idx = hist2d_chisq_corr > 0
    # img.scatter_pts(0,2,xh[idx], yh[idx], hist2d_diff_corr[idx],color_map="jet")
    # img.scatter_pts(1,2,xh[idx], yh[idx], hist2d_chisq_corr[idx],color_map="jet")
    #
    # img.save_img("./result/hist2d/hist2d_%d_%d.png"%(tag_i,tag_j))
    # img.close_img()


    img = Image_Plot(xpad=0.25, ypad=0.25)
    img.subplots(2,3)

    img.axs[0][0].scatter(mg_pdf_bin_label[mg_pdf_bin_num2:], hist1d_diff_ori)
    img.axs[1][0].scatter(mg_pdf_bin_label[mg_pdf_bin_num2:], hist1d_chisq_ori)

    img.axs[0][1].scatter(mg_pdf_bin_label[mg_pdf_bin_num2:], hist1d_diff_ori_sym)
    img.axs[1][1].scatter(mg_pdf_bin_label[mg_pdf_bin_num2:], hist1d_chisq_ori_sym)

    img.axs[0][2].scatter(mg_pdf_bin_label[mg_pdf_bin_num2:], hist1d_diff_corr)
    img.axs[1][2].scatter(mg_pdf_bin_label[mg_pdf_bin_num2:], hist1d_chisq_corr)

        #
        # for j in range(3):
        #     if i == 0:
        #         img.axs[i][j].set_title("before PDF_SYM")
        #     else:
        #         img.axs[i][j].set_title("after PDF_SYM")
        #     img.set_label(i,j, 1, "bin label")
        #     img.axs[i][j].legend()
        #     img.axis_sci_ticklabel(i, j, 0)
        # img.set_label(i, 1, 0, "$n_{right} - n_{left}$")
        # img.set_label(i, 2, 0, "$\chi^2$")

    img.save_img("./result/hist1d/hist1d_%s_%d_%d.png"%(data_type,tag_i,tag_j))
    img.close_img()
    #
comm.Barrier()
if rank == 0:
    numpy.savez("./result/cache_%s.npz"%data_type, final_result, delta_sigma)
    img = Image_Plot(xpad=0.25,ypad=0.25)
    img.subplots(2,2)

    img.scatter_pts(0,0,scale_gridx, scale_gridy, final_result[:element_num].flatten(),color_map="jet")
    img.scatter_pts(0,1,scale_gridx, scale_gridy, final_result[element_num: int(2*element_num)].flatten(),color_map="jet")

    img.scatter_pts(1,0,scale_gridx, scale_gridy, final_result[int(2*element_num): int(3*element_num)].flatten(),color_map="jet")
    img.axs[1][0].set_title("True $\Delta\Sigma$=%.4f"%delta_sigma)
    img.scatter_pts(1,1,scale_gridx, scale_gridy, final_result[int(3*element_num): int(4*element_num)].flatten(),color_map="jet")

    img.save_img("./result/corr_result_%s.png"%data_type)

    img.close_img()