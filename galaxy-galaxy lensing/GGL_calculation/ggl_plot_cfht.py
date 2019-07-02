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


argc = len(argv)
area_num = argc - 3

fore_source = argv[1]
cmd = argv[2]

total_path = "/mnt/perc/hklee/CFHT/gg_lensing/"

h5f = h5py.File(total_path + "result/%s/cfht/w_%d/radius_bin.hdf5"%(fore_source, 1), "r")
radius_bin = h5f["/radius_bin"].value[:,0]
h5f.close()
radius_num = radius_bin.shape[0]-1

# radius_num = 11
# radius_bin = tool_box.set_bin_log(0.1, 15, radius_num+1)


if area_num > 1:
    result_path = total_path + "result/%s/cfht/result/cfht_%s_result_total.hdf5"%(fore_source,fore_source)
    dens_pic_path = total_path + "result/%s/cfht/result/cfht_%s_total"%(fore_source,fore_source)
    dens_r_pic_path = total_path + "result/%s/cfht/result/cfht_%s_total_sgimaxr"%(fore_source,fore_source)
else:
    result_path = total_path + "result/%s/cfht/result/cfht_%s_result_w_%d.hdf5"%(fore_source, fore_source, int(argv[3]))
    dens_pic_path = total_path + "result/%s/cfht/result/cfht_%s_w_%d"%(fore_source, fore_source, int(argv[3]))
    dens_r_pic_path = total_path + "result/%s/cfht/result/cfht_%s_w_%d_sigmaxr"%(fore_source, fore_source, int(argv[3]))

ylabels = ["$\gamma$", "$\Delta\Sigma \; [\\rm{M_{\odot}} \cdot \\rm{pc^{-2}}]$"]
ylabels_r = "$\\rm{R}\Delta\Sigma \; [\\rm{10^6\ M_{\odot}} \cdot \\rm{pc^{-2}}]$"
xlabel = "$\\rm{R} \; [\\rm{Mpc} \cdot \\rm{h^{-1}}]$"

coeff = 388.2833518

gt_lb, gx_lb = 0, 2
sigt_lb, sigx_lb = 4, 6
sigtxr_lb = 8
r_lb = 10

if cmd == "calculate":
    result = numpy.zeros((11, radius_num))
    for ir in range(radius_num):

        stack_count = 0
        for ia in range(3, argc):
            h5f = h5py.File(total_path + "result/%s/cfht/w_%d/radius_%d.hdf5"%(fore_source, int(argv[ia]), ir), "r")
            temp = h5f["/pair_data"].value
            h5f.close()

            if temp[0,8] > -1:
                if stack_count == 0:
                    data = temp
                else:
                    data = numpy.row_stack((data, temp))
                stack_count += 1
        if stack_count > 0:
            pair_num = data.shape[0]
            e1 = data[:, 0]
            # additive bias c2
            e2 = -data[:, 1]# + data[:, 5]
            cos_2phi = data[:, 2]
            sin_2phi = data[:, 3]
            m_bias = data[:, 4]
            weight_measure = data[:, 6]
            crit = data[:, 7]
            # radius from the center
            dist = data[:,8]

            weight = weight_measure/crit**2
            weight_sum = weight.sum()

            e_t = e1*cos_2phi - e2*sin_2phi
            e_x = e1*sin_2phi + e2*cos_2phi

            # two kinds of correction
            corr_m = 1 + numpy.sum(weight * m_bias)/weight_sum
            corr_m_ = 1 + m_bias.mean()

            delta_crit_et = e_t*crit*coeff
            delta_crit_ex = e_x*crit*coeff
            delta_sigma_t = numpy.sum(weight*delta_crit_et)/weight_sum/corr_m_
            delta_sigma_x = numpy.sum(weight*delta_crit_ex)/weight_sum/corr_m_
            gamma_t = numpy.sum(weight*e_t)/weight_sum/corr_m_
            gamma_x = numpy.sum(weight*e_x)/weight_sum/corr_m_

            r_mean = dist.mean()

            result[gt_lb, ir] = gamma_t
            result[gx_lb, ir] = gamma_x
            result[sigt_lb, ir] = delta_sigma_t
            result[sigx_lb, ir] = delta_sigma_x
            result[sigtxr_lb, ir] = delta_sigma_t*r_mean

            result[gt_lb+1, ir] = e_t.std()/numpy.sqrt(pair_num)
            result[gx_lb+1, ir] = e_x.std()/numpy.sqrt(pair_num)
            result[sigt_lb+1, ir] = delta_crit_et.std()/numpy.sqrt(pair_num)
            result[sigx_lb+1, ir] = delta_crit_ex.std()/numpy.sqrt(pair_num)
            result[sigtxr_lb+1, ir] = delta_crit_et.std()/numpy.sqrt(pair_num)*r_mean

            result[r_lb, ir] = r_mean

            print("[%.5f, %.5f], %d galaxy pairs at radius %f (%f)"%(radius_bin[ir],radius_bin[ir+1], pair_num,
                                                                     r_mean,(radius_bin[ir]+radius_bin[ir+1])/2))
            # tag = numpy.arange(0, pair_num)
            #
            # for jnf in range(50):
            #     ch = numpy.random.choice(tag,pair_num)
            #
            #     weight_ch = weight[tag]
            #     weight_ch_sum = weight_ch.sum()
            #     crit_ch = crit[tag]
            #     e_t_ch = e_t[tag]
            #     e_x_ch = e_x[tag]
            #     n_bias_ch = m_bias[tag]
            #
            #     delta_sigma_t = numpy.sum(weight * e_t * crit) / weight_sum / corr_m * 388.2833518
            #     delta_sigma_x = numpy.sum(weight * e_x * crit) / weight_sum / corr_m * 388.2833518
            #     gamma_t = numpy.sum(weight * e_t) / weight_sum / corr_m
            #     gamma_x = numpy.sum(weight * e_x) / weight_sum / corr_m
        else:
            print("Skip [%.5f, %.5f], 0 galaxy pairs"%(radius_bin[ir],radius_bin[ir+1]))
    h5f = h5py.File(result_path,"w")
    h5f["/data"] = result
    h5f.close()

    ylims = (0.01, 140)

    img = Image_Plot()
    img.set_style()
    img.subplots(1, 1)
    # img.axs[0][0].errorbar(result[r_lb], result[gt_lb], result[gt_lb + 1], c="C1", capsize=4, label="T", marker="s")
    # img.axs[0][0].errorbar(result[r_lb], result[gx_lb], result[gx_lb + 1], c="C2", capsize=4, label="X", marker="s")

    img.axs[0][0].errorbar(result[r_lb], result[sigt_lb], result[sigt_lb + 1], c="C1", capsize=4, label="T", marker="s")
    img.axs[0][0].errorbar(result[r_lb], result[sigx_lb], result[sigx_lb + 1], c="C2", capsize=4, label="X", marker="s")

    # plot the line of "W1" extracted from "Lensing is low"
    if area_num == 1 and int(argv[3]) == 1 and fore_source == "cmass":
        w1_cfht_path = "./lensing_low/data.dat"
        if os.path.exists(w1_cfht_path):
            w1_data_cfht = numpy.loadtxt(w1_cfht_path)
            img.axs[0][0].errorbar(w1_data_cfht[:, 0], w1_data_cfht[:, 1], w1_data_cfht[:, 2], c="C4", capsize=4, label="w1, Lensing_low", marker="s")

    img.set_label(0, 0, 0, ylabels[1])
    img.set_label(0, 0, 1, xlabel)

    img.axs[0][0].set_yscale("log")
    # img.axs[0][0].set_ylim(ylims)
    img.axs[0][0].set_xscale("log")
    xs = img.axs[0][0].set_xlim()
    img.axs[0][0].plot([xs[0], xs[1]], [0, 0], linestyle="--", linewidth=1, c="grey")
    img.set_legend(0,0,loc="upper right")

    for j in range(10):
        img.axs[0][0].plot([xs[0], xs[1]], [j, j], linewidth=0.5, c="grey", alpha=0.5)
        img.axs[0][0].plot([xs[0], xs[1]], [10 + 10*j, 10 + 10*j], linewidth=0.5,c="grey", alpha=0.5)
        img.axs[0][0].plot([xs[0], xs[1]], [100 + 100*j, 100 + 100*j], linewidth=0.5,c="grey", alpha=0.5)

    img.axs[0][0].set_xlim(xs[0], xs[1])

    img.save_img(dens_pic_path + ".png")
    img.set_style_default()
    img.close_img()

    # plot R x \Delta\Sigma
    img = Image_Plot()
    img.set_style()
    img.subplots(1,1)
    img.set_label(0, 0, 0, ylabels_r)
    img.set_label(0, 0, 1, xlabel)
    img.axs[0][0].errorbar(result[r_lb], result[sigtxr_lb], result[sigtxr_lb + 1], c="C1", capsize=4, label="X", marker="s")
    img.axs[0][0].set_xscale("log")
    img.save_img(dens_r_pic_path + ".png")
    img.set_style_default()
    img.close_img()

if cmd == "plot":

    h5f = h5py.File(result_path,"r")
    result = h5f["/data"].value
    h5f.close()

    ylims = (0.01, 140)

    img = Image_Plot()
    img.set_style()
    img.subplots(1, 1)
    # img.axs[0][0].errorbar(result[r_lb], result[gt_lb], result[gt_lb + 1], c="C1", capsize=4, label="T", marker="s")
    # img.axs[0][0].errorbar(result[r_lb], result[gx_lb], result[gx_lb + 1], c="C2", capsize=4, label="X", marker="s")

    img.axs[0][0].errorbar(result[r_lb], result[sigt_lb], result[sigt_lb + 1], c="C1", capsize=4, label="T", marker="s")
    img.axs[0][0].errorbar(result[r_lb], result[sigx_lb], result[sigx_lb + 1], c="C2", capsize=4, label="X", marker="s")

    # plot the line of "W1" extracted from "Lensing is low"
    if area_num == 1 and int(argv[3]) == 1 and fore_source == "cmass":
        w1_cfht_path = "./lensing_low/data.dat"
        if os.path.exists(w1_cfht_path):
            w1_data_cfht = numpy.loadtxt(w1_cfht_path)
            img.axs[0][0].errorbar(w1_data_cfht[:, 0], w1_data_cfht[:, 1], w1_data_cfht[:, 2], c="C4", capsize=4, label="w1, Lensing_low", marker="s")

    img.set_label(0, 0, 0, ylabels[1])
    img.set_label(0, 0, 1, xlabel)

    img.axs[0][0].set_yscale("log")
    # img.axs[0][0].set_ylim(ylims)
    img.axs[0][0].set_xscale("log")
    xs = img.axs[0][0].set_xlim()
    img.axs[0][0].plot([xs[0], xs[1]], [0, 0], linestyle="--", linewidth=1, c="grey")
    img.set_legend(0,0,loc="upper right")

    for j in range(10):
        img.axs[0][0].plot([xs[0], xs[1]], [j, j], linewidth=0.7, c="grey", alpha=0.6)
        img.axs[0][0].plot([xs[0], xs[1]], [10 + 10*j, 10 + 10*j], linewidth=0.7,c="grey", alpha=0.6)
        img.axs[0][0].plot([xs[0], xs[1]], [100 + 100*j, 100 + 100*j], linewidth=0.7,c="grey", alpha=0.6)

    img.axs[0][0].set_xlim(xs[0], xs[1])

    img.save_img(dens_pic_path + ".png")
    img.set_style_default()
    img.close_img()

    # plot R x \Delta\Sigma
    img = Image_Plot()
    img.set_style()
    img.subplots(1,1)
    img.set_label(0, 0, 0, ylabels_r)
    img.set_label(0, 0, 1, xlabel)
    img.axs[0][0].errorbar(result[r_lb], result[sigtxr_lb], result[sigtxr_lb + 1], c="C1", capsize=4, label="X", marker="s")
    img.axs[0][0].set_xscale("log")
    img.save_img(dens_r_pic_path + ".png")
    img.set_style_default()
    img.close_img()

print("Images are saved in %s"%dens_pic_path)

