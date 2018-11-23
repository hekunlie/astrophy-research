from subprocess import Popen

sources = ["dimmer", "dimmerm3"]

for source in sources:

    cmd = "mpirun -np 60 python /home/hkli/work/selection_bias/CFHT_simu/selection_bias_images_h.py %s"%source
    a = Popen(cmd, shell=True)
    a.wait()

    # cmd = "mpiexec -n 14 /home/hkli/work/c/measure_dimmer"
    # a = Popen(cmd, shell=True)
    # a.wait()

    # if "pts" in source:
    #     max_radius = 8
    # else:
    #     max_radius = 5.5
    #
    # filter_names = ["sex2_2", "sex3_2", "sex4_2"]
    # gauss_filters = ["gauss_2.0_5x5", "gauss_3.0_5x5", "gauss_4.0_7x7"]
    #
    # for ii, filter_name in enumerate(filter_names):
    #     # change the filter
    #     with open("./default.sex_ori", "r") as f:
    #         contents = f.readlines()
    #     contents[23] = "FILTER_NAME      %s.conv   # name of the file containing the filter"%gauss_filters[ii]
    #     with open("./default.sex", "w") as f:
    #         f.writelines(contents)
    #
    #     # run the estimation code
    #     cmd = "mpirun -np 28 python snr_est.py snr %s %s %.1f"%(source, filter_name, max_radius)
    #     a = Popen(cmd, shell=True)
    #     a.wait()
    #
    #     # add the data in .cat files to the final catalogs
    #     cmd = "mpirun -np 14 python snr_est.py add %s %s %.1f"%(source, filter_name, max_radius)
    #     a = Popen(cmd, shell=True)
    #     a.wait()
    #
    #     # check
    #     cmd = "mpirun -np 14 python snr_est.py check %s %s %.1f"%(source, filter_name, max_radius)
    #     a = Popen(cmd, shell=True)
    #     a.wait()

# exit()
# file_num = 1
# resolution_factor = [0.7]
#
# cuts = ['sex_snr', 'mag_auto', "flux2"]
# filter_name = ["sex2_2", "sex3_2","sex4_2"]
# sig = ["2sig", "2sig", "2sig"]
#
# # cuts = ['sex_snr', 'mag_auto', "flux2"]
# # filter_name = ["sex2_1.5", "sex2_2", "sex3_1.5", "sex3_2", "sex4_1.5", "sex4_2"]
# # sig = ["1.5sig", "2sig", "1.5sig", "2sig", "1.5sig", "2sig"]
# for s, source in enumerate(sources):
#     for f in range(len(filter_name)):
#         for i in range(len(cuts)):
#             if i < 2 or "sex2_" in filter_name[f]:
#                 cmd = "mpirun -np 14 python /home/hkli/work/selection_bias/sym_mc_plot/" \
#                       "sym_mc_plot.py %s %d %s %s %s %.2f"%\
#                       (cuts[i], file_num, filter_name[f], source, sig[f], resolution_factor[s])
#                 a = Popen(cmd, shell=True)
#                 a.wait()

