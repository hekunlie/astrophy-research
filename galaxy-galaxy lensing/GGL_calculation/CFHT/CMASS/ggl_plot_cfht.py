import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import matplotlib.pyplot as plt
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from plot_tool import Image_Plot
import tool_box
import h5py


# plot the GGL signal
# "python .. cmass plot 1" or "python .. cmass calculate 1 3 4

argc = len(argv)
area_num = argc - 3
# foreground
fore_source = argv[1]
# "calculate" or "plot"
cmd = argv[2]

parent_path = "/mnt/perc/hklee/CFHT/gg_lensing/"
parent_result_path = parent_path + "result/%s/cfht/"%fore_source

lens_low_path = '/home/hklee/work/CFHT/gg_lensing/lensing_low/data.dat'

h5f = h5py.File(parent_result_path + "w_%s/radius_0.hdf5"%argv[3], "r")
radius_bin = h5f["/radius_bin"].value[:,0]
h5f.close()
radius_num = radius_bin.shape[0]-1

# radius_num = 11
# radius_bin = tool_box.set_bin_log(0.1, 15, radius_num+1)


if area_num > 1:
    result_path = parent_result_path + "result/%s_result_total.hdf5"%fore_source
    result_path_npz = parent_result_path + "result/%s_result_total.npz"%fore_source
    dens_pic_path = parent_result_path + "result/%s_total"%fore_source
    dens_r_pic_path = parent_result_path + "result/%s_total_sgimaxr"%fore_source
else:
    result_path = parent_result_path + "result/%s_result_w_%d.hdf5"%(fore_source, int(argv[3]))
    result_path_npz = parent_result_path + "result/%s_result_w_%d.npz"%(fore_source, int(argv[3]))
    dens_pic_path = parent_result_path + "result/%s_w_%d"%(fore_source, int(argv[3]))
    dens_r_pic_path = parent_result_path + "result/%s_w_%d_sigmaxr"%(fore_source, int(argv[3]))

ylabels = ["$\gamma$", "$\Delta\Sigma \; [\\rm{h \cdot M_{\odot}} \cdot \\rm{pc^{-2}}]$"]
ylabels_r = "$\\rm{R}\Delta\Sigma \; [\\rm{10^6\ M_{\odot}} \cdot \\rm{pc^{-2}}]$"
xlabel = "$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$"

coeff_integ = 554.6821355281792
coeff_dist = 1662895.2081868195

# the catalog contains
# e_t, e_x, m, c, weight, crit_density,
# integrate part of crit_density, transverse distance and redshift of background

# result label
crit_t_lb, crit_t_sig_lb = 0, 1
crit_x_lb, crit_x_sig_lb = 2, 3
trans_dist_lb = 4

if cmd == "calculate":
    result = numpy.zeros((5, radius_num))
    for ir in range(radius_num):

        stack_count = 0
        for ia in range(3, argc):
            h5f = h5py.File(parent_result_path + "w_%d/radius_%d.hdf5"%(int(argv[ia]), ir), "r")
            try:
                sub_num = h5f["/subset_num"].value[0,0]

                for i_sub in range(sub_num):
                    temp = h5f["/pair_data_%d"%i_sub].value

                    if temp.shape[0] > 0:
                        if stack_count == 0:
                            data = temp
                        else:
                            data = numpy.row_stack((data, temp))
                        stack_count += 1
            except:
                pass
            h5f.close()

        if stack_count > 0:

            pair_num = data.shape[0]

            e_t = data[:, 0]
            e_x = data[:, 1]
            # m bias
            m_bias = data[:, 2]
            # weight
            weight_measure = data[:, 3]
            crit_integ = data[:, 4]
            crit = data[:, 5]
            # radius from the center
            dist = data[:, 6]
            redshif = data[:, 7]

            # weight = weight_measure/crit_integ**2
            weight = weight_measure/crit**2
            weight_sum = weight.sum()
            # weight_bias = weight * m_bias
            # weight_sum = tool_box.accurate_sum(weight, 10000)

            # correction
            corr_m = 1 + numpy.sum(weight * m_bias)/weight_sum

            # delta_crit_et = e_t*crit_integ*coeff_integ*weight
            # delta_crit_ex = e_x*crit_integ*coeff_integ*weight

            delta_crit_et = e_t*crit*coeff_dist*weight
            delta_crit_ex = e_x*crit*coeff_dist*weight

            delta_sigma_t = numpy.sum(delta_crit_et)/weight_sum/corr_m
            err_t = numpy.sum(e_t.std()*crit*coeff_dist*weight)/weight_sum/corr_m/numpy.sqrt(pair_num)

            delta_sigma_x = numpy.sum(delta_crit_ex)/weight_sum/corr_m
            err_x = numpy.sum(e_x.std()*crit*coeff_dist*weight)/weight_sum/corr_m/numpy.sqrt(pair_num)

            # delta_sigma_t = tool_box.accurate_sum(delta_crit_et, 10000)/weight_sum/corr_m
            # delta_sigma_x = tool_box.accurate_sum(delta_crit_ex, 10000)/weight_sum/corr_m

            # gamma_t = numpy.sum(weight*e_t)/weight_sum/corr_m
            # gamma_x = numpy.sum(weight*e_x)/weight_sum/corr_m

            r_mean = dist.mean()
            # r_mean = tool_box.accurate_sum(dist, 1000)/dist.shape[0]

            result[crit_t_lb, ir] = delta_sigma_t
            result[crit_t_sig_lb, ir] = err_t
            result[crit_x_lb, ir] = delta_sigma_x
            result[crit_x_sig_lb, ir] = err_x
            result[trans_dist_lb, ir] = r_mean

            print("[%.5f, %.5f], %d galaxy pairs at radius %f (%f). ESD: %.3f (%.3f)"%(
                radius_bin[ir],radius_bin[ir+1], pair_num, r_mean,(radius_bin[ir]+radius_bin[ir+1])/2,result[crit_t_lb, ir],result[crit_t_sig_lb, ir]))
        else:
            print("Skip [%.5f, %.5f], 0 galaxy pairs"%(radius_bin[ir],radius_bin[ir+1]))
    h5f = h5py.File(result_path,"w")
    h5f["/data"] = result
    h5f.close()
    numpy.savez(result_path_npz, result)

    img = Image_Plot()
    img.set_style()
    img.subplots(1, 1)
    # img.axs[0][0].errorbar(result[r_lb], result[gt_lb], result[gt_lb + 1], c="C1", capsize=4, label="T", marker="s")
    # img.axs[0][0].errorbar(result[r_lb], result[gx_lb], result[gx_lb + 1], c="C2", capsize=4, label="X", marker="s")

    img.axs[0][0].errorbar(result[trans_dist_lb], result[crit_t_lb], result[crit_t_sig_lb],
                            c="C1", marker="s", capsize=4, mfc="none",fmt=" ",label="T")
    img.axs[0][0].errorbar(result[trans_dist_lb], result[crit_x_lb], result[crit_x_sig_lb + 1],
                           c="C2", marker="s", capsize=4,mfc="none", fmt=" ",label="X")

    y_max = img.axs[0][0].set_ylim()[1]
    ylims = (0.01, 200)
    # plot the line of "W1" extracted from "Lensing is low"
    if area_num == 1 and fore_source == "cmass":
        if os.path.exists(lens_low_path):
            w1_data_cfht = numpy.loadtxt(lens_low_path)
            img.axs[0][0].errorbar(w1_data_cfht[:, 0], w1_data_cfht[:, 1], w1_data_cfht[:, 2], marker="s",
                                   c="C4", capsize=4, mfc="none",fmt=" ", label="w1, Lensing_low")

    img.set_label(0, 0, 0, ylabels[1])
    img.set_label(0, 0, 1, xlabel)

    img.axs[0][0].set_yscale("log")
    img.axs[0][0].set_ylim(ylims)
    img.axs[0][0].set_xscale("log")
    xs = img.axs[0][0].set_xlim()
    # img.axs[0][0].plot([xs[0], xs[1]], [0, 0], linestyle="--", linewidth=1, c="grey")
    img.set_legend(0,0,loc="upper right")

    # for j in range(10):
    #     img.axs[0][0].plot([xs[0], xs[1]], [j, j], linewidth=0.5, c="grey", alpha=0.5)
    #     img.axs[0][0].plot([xs[0], xs[1]], [10 + 10*j, 10 + 10*j], linewidth=0.5,c="grey", alpha=0.5)
    #     img.axs[0][0].plot([xs[0], xs[1]], [100 + 100*j, 100 + 100*j], linewidth=0.5,c="grey", alpha=0.5)
    #
    # img.axs[0][0].set_xlim(xs[0], xs[1])

    img.save_img(dens_pic_path + ".png")
    img.set_style_default()
    img.close_img()

    # # plot R x \Delta\Sigma
    # img = Image_Plot()
    # img.set_style()
    # img.subplots(1,1)
    # img.set_label(0, 0, 0, ylabels_r)
    # img.set_label(0, 0, 1, xlabel)
    # img.axs[0][0].errorbar(result[trans_dist_lb], result[sigtxr_lb], result[sigtxr_lb + 1], c="C1", capsize=4, label="X", marker="s")
    # img.axs[0][0].set_xscale("log")
    # img.save_img(dens_r_pic_path + ".png")
    # img.set_style_default()
    # img.close_img()

if cmd == "plot":

    h5f = h5py.File(result_path,"r")
    result = h5f["/data"].value
    h5f.close()

    img = Image_Plot()
    img.set_style()
    img.subplots(1, 1)
    # img.axs[0][0].errorbar(result[r_lb], result[gt_lb], result[gt_lb + 1], c="C1", capsize=4, label="T", marker="s")
    # img.axs[0][0].errorbar(result[r_lb], result[gx_lb], result[gx_lb + 1], c="C2", capsize=4, label="X", marker="s")

    img.axs[0][0].errorbar(result[trans_dist_lb], result[crit_t_lb], result[crit_t_sig_lb], c="C1", mfc="none", marker="s",
                           capsize=4,fmt=" ", label="T")
    img.axs[0][0].errorbar(result[trans_dist_lb], result[crit_x_lb], result[crit_x_sig_lb + 1], c="C2", mfc="none", marker="s",
                           capsize=4,fmt=" ", label="X")

    y_max = img.axs[0][0].set_ylim()[1]
    ylims = (0.01, 200)

    # plot the line of "W1" extracted from "Lensing is low"
    if area_num == 1 and int(argv[3]) == 1 and fore_source == "cmass":
        if os.path.exists(lens_low_path):
            w1_data_cfht = numpy.loadtxt(lens_low_path)
            img.axs[0][0].errorbar(w1_data_cfht[:, 0], w1_data_cfht[:, 1], w1_data_cfht[:, 2], c="C4", mfc="none",
                                   capsize=4, marker="s", label="w1, Lensing_low",fmt=" ")

    img.set_label(0, 0, 0, ylabels[1])
    img.set_label(0, 0, 1, xlabel)

    img.axs[0][0].set_yscale("log")
    img.axs[0][0].set_ylim(ylims)
    img.axs[0][0].set_xscale("log")
    xs = img.axs[0][0].set_xlim()
    # img.axs[0][0].plot([xs[0], xs[1]], [0, 0], linestyle="--", linewidth=1, c="grey")
    img.set_legend(0,0,loc="upper right")

    # for j in range(10):
    #     img.axs[0][0].plot([xs[0], xs[1]], [j, j], linewidth=0.7, c="grey", alpha=0.6)
    #     img.axs[0][0].plot([xs[0], xs[1]], [10 + 10*j, 10 + 10*j], linewidth=0.7,c="grey", alpha=0.6)
    #     img.axs[0][0].plot([xs[0], xs[1]], [100 + 100*j, 100 + 100*j], linewidth=0.7,c="grey", alpha=0.6)
    #
    # img.axs[0][0].set_xlim(xs[0], xs[1])

    img.save_img(dens_pic_path + ".png")
    img.set_style_default()
    img.close_img()

    # # plot R x \Delta\Sigma
    # img = Image_Plot()
    # img.set_style()
    # img.subplots(1,1)
    # img.set_label(0, 0, 0, ylabels_r)
    # img.set_label(0, 0, 1, xlabel)
    # img.axs[0][0].errorbar(result[r_lb], result[sigtxr_lb], result[sigtxr_lb + 1], c="C1", capsize=4, label="X", marker="s")
    # img.axs[0][0].set_xscale("log")
    # img.save_img(dens_r_pic_path + ".png")
    # img.set_style_default()
    # img.close_img()

print("Images are saved in %s"%dens_pic_path)

