import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/fourier_quad/'%my_home)
path.append('E:/Github/astrophy-research/my_lib/')
import tool_box
import shutil


source, cmd = argv[1], argv[2]

envs_path = "%s/work/envs/envs.dat"%my_home
get_contents = [['selection_bias', "%s_path"%source, '1'],['selection_bias', "%s_path_result"%source, '1'],
                ['selection_bias', "%s_path_para"%source, '1'],['selection_bias', "%s_path_log"%source, '1']]
path_items = tool_box.config(envs_path,['get','get','get','get'], get_contents)
total_path, result_path, para_path, log_path = path_items


sex_filters = ["sex2_2", "sex3_2", "sex4_2", "sex2_1.5", "sex3_1.5", "sex4_1.5"]
cut_nm = ["mag_auto", "sex_snr", "flux2"]
if cmd == "cut":
    cut_path = result_path + "cuts/"
    if os.path.exists(cut_path):
        shutil.rmtree(cut_path)
    for sub_filter in sex_filters:
        for sub_nm in cut_nm:
            sub_path = cut_path + "sym/%s/%s/"%(sub_filter, sub_nm)
            os.makedirs(sub_path)
            print("Build: %s"%sub_path)


if cmd == "sex":
    for sub_sex in sex_filters:
        sex_path = result_path + "data/%s/"%sub_sex
        if os.path.exists(sex_path):
            shutil.rmtree(sex_path)
        os.makedirs(sex_path + "cat/")
        print("Build: %s/cat/" % sex_path)

if cmd == "all":
    for i in range(14):
        img_path = total_path + "%d/"
        os.makedirs(img_path)
    os.makedirs(para_path + "logs/")
    os.makedirs(para_path + "pic/")
    os.mkdir(total_path + "logs/")
    os.makedirs(result_path + "data/data_2sig/")
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