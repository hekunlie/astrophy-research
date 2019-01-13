import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
path.append('E:/Github/astrophy-research/my_lib/')
import tool_box
import shutil


sect, source, cmd = argv[1], argv[2], argv[3]

envs_path = "%s/work/envs/envs.dat"%my_home
get_contents = [['%s'%sect, "%s_path"% source, '1'],['%s'%sect, "%s_path_result"%source, '1'],
                ['%s'%sect, "%s_path_para"%source, '1'],['%s'%sect, "%s_path_log"%source, '1']]
path_items = tool_box.config(envs_path,['get','get','get','get'], get_contents)
total_path, result_path, para_path, log_path = path_items


sex_filters = ["sex2_2", "sex3_2", "sex4_2", "sex2_1.5", "sex3_1.5", "sex4_1.5"]
cut_nm = ["mag_auto", "sex_snr", "flux2", "flux_alt", "snr_auto", "flux", "snr",
          "flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5"]

if cmd == "cut":
    cut_path = result_path + "cuts/"
    if os.path.exists(cut_path):
        shutil.rmtree(cut_path)
    for sub_filter in sex_filters:
        for sub_nm in cut_nm:
            sub_path = cut_path + "sym/%s/%s/"%(sub_filter, sub_nm)
            if "sex2_" in sub_filter:
                os.makedirs(sub_path)
                print("Build: %s"%sub_path)
            else:
                if "mag" in sub_nm or "sex" in sub_nm or "snr_" in sub_nm:
                    os.makedirs(sub_path)
                    print("Build: %s" % sub_path)


if cmd == "sex":
    for sub_sex in sex_filters:
        sex_path = result_path + "data/%s/"%sub_sex
        if os.path.exists(sex_path):
            shutil.rmtree(sex_path)
        os.makedirs(sex_path + "cat/")
        print("Build: %scat/" % sex_path)

if cmd == "data":
    os.makedirs(result_path + "data/data_2.0sig/")
    print("Build: %scat/" % (result_path + "data/data_2.0sig/"))
    os.makedirs(result_path + "data/data_1.5sig/")
    print("Build: %scat/" % (result_path + "data/data_1.5sig/"))
    os.makedirs(result_path + "data/check/")
    print("Build: %scat/" % (result_path + "data/check/"))
    for sub_sex in sex_filters:
        sex_path = result_path + "data/%s/"%sub_sex
        if os.path.exists(sex_path):
            shutil.rmtree(sex_path)
        os.makedirs(sex_path + "cat/")
        print("Build: %scat/" % sex_path)

if cmd == "all":
    for i in range(14):
        img_path = total_path + "%d/"%i
        os.makedirs(img_path)
    os.makedirs(para_path + "logs/")
    os.makedirs(para_path + "pic/")
    os.mkdir(total_path + "logs/")
    os.makedirs(result_path + "data/data_2.0sig/")
    os.makedirs(result_path + "data/data_1.5sig/")
    os.makedirs(result_path + "data/check/")
    os.makedirs(result_path + "pic/")

    cut_path = result_path + "cuts/"
    for sub_sex in sex_filters:
        sex_path = result_path + "data/%s/"%sub_sex
        os.makedirs(sex_path + "cat/")
        sex_check = result_path + "data/check/%s/"%sub_sex
        os.makedirs(sex_check)
        for sub_nm in cut_nm:
            sub_path = cut_path + "sym/%s/%s/"%(sub_sex, sub_nm)
            os.makedirs(sub_path)