import os
from sys import path, argv
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
path.append('%s/work/mylib/'%my_home)
path.append('E:/Github/astrophy-research/mylib/')
import tool_box
import shutil


total_path = argv[1]
cmd = argv[2]


result_path = total_path + "/result/"
para_path = total_path + "/parameters/"
log_path = total_path + "/logs/"


sex_filters = ["sex2_4", "sex3_4", "sex4_4","sex2_2", "sex3_2", "sex4_2", "sex2_1.5", "sex3_1.5", "sex4_1.5"]
cut_nm = ["mag_auto", "snr_sex", "snr_auto",  "flux2_ex1", "flux2_ex2", "flux2_ex3", "flux2_ex4", "flux2_ex5","rfactor"]

if cmd == "cut":
    cut_path = result_path + "cuts/"
    if os.path.exists(cut_path):
        shutil.rmtree(cut_path)
    for sub_filter in sex_filters:
        for sub_nm in cut_nm:
            sub_path = cut_path + "sym/%s/%s/"%(sub_filter, sub_nm)
            os.makedirs(sub_path)
            print("Build: %s" % sub_path)
            # if "sex2_" in sub_filter:
            #     os.makedirs(sub_path)
            #     print("Build: %s"%sub_path)
            # else:
            #     if "mag" in sub_nm or "sex" in sub_nm or "snr_" in sub_nm:
            #         os.makedirs(sub_path)
            #         print("Build: %s" % sub_path)


if cmd == "sex":
    for sub_sex in sex_filters:
        sex_path = result_path + "data/%s/"%sub_sex
        if not os.path.exists(sex_path):
            # shutil.rmtree(sex_path)
            os.makedirs(sex_path + "cat/")
        print("Build: %scat/" % sex_path)

if cmd == "data":
    data_sigs = ["data/data_2.0sig/","data/data_1.5sig/","data/data_4.0sig/"]
    for sigs in data_sigs:
        data_path = result_path + sigs
        if not os.path.exists(data_path):
            #shutil.rmtree(data_path)
            os.makedirs(data_path)
            print("Build: %s/" %data_path)
        else:
            print(data_path," exists")
    for sub_sex in sex_filters:
        sex_path = result_path + "data/%s/"%sub_sex
        if not os.path.exists(sex_path):
            # shutil.rmtree(sex_path)
            os.makedirs(sex_path + "cat/")
        print("Build: %scat/" % sex_path)

if cmd == "all":
    for i in range(14):
        img_path = total_path + "%d/"%i
        if not os.path.exists(img_path):
            os.makedirs(img_path)

    if not os.path.exists(para_path + "logs/"):
        os.makedirs(para_path + "logs/")

    if not os.path.exists(para_path + "pic/"):
        os.makedirs(para_path + "pic/")

    if not os.path.exists(total_path + "logs/"):
        os.makedirs(total_path + "logs/")

    if not os.path.exists(result_path + "data/data_2.0sig/"):
        os.makedirs(result_path + "data/data_2.0sig/")
    if not os.path.exists(result_path + "data/data_1.5sig/"):
        os.makedirs(result_path + "data/data_1.5sig/")
    if not os.path.exists(result_path + "data/data_4.0sig/"):
        os.makedirs(result_path + "data/data_4.0sig/")

    if not os.path.exists(result_path + "data/check/"):
        os.makedirs(result_path + "data/check/")

    if not os.path.exists(result_path + "pic/"):
        os.makedirs(result_path + "pic/")

    cut_path = result_path + "cuts/"
    for sub_sex in sex_filters:
        sex_path = result_path + "data/%s/"%sub_sex
        if not os.path.exists(sex_path + "cat/"):
            os.makedirs(sex_path + "cat/")

        sex_check = result_path + "data/check/%s/" % sub_sex
        if not os.path.exists(sex_check):
            os.makedirs(sex_check)

        for sub_nm in cut_nm:
            sub_path = cut_path + "sym/%s/%s/"%(sub_sex, sub_nm)
            if not os.path.exists(sub_path):
                os.makedirs(sub_path)