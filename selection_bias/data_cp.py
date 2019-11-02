from subprocess import Popen
import os

sources = ["galsim_bright", "galsim_dimmer", "pts_bright", "pts_dimmer"]
sex_nms = ["sex2_1.5", "sex2_2", "sex2_4", "sex4_1.5", "sex4_2", "sex4_4"]
data_nms = []
for i in range(10):
    data_nms.append("Rfacotr_%d.hdf5"%i)
    data_nms.append("sex_%d.hdf5"%i)
    data_nms.append("mask_%d.hdf5"%i)

# for src_nm in sources:
#     for src_sex in sex_nms:
#         dst_path = "/mnt/perc/hklee/selection_bias/%s/result/data/%s"%(src_nm, src_sex)
#         if not os.path.exists(dst_path):
#             os.mkdir(dst_path)
#         for dd in data_nms:
#             cmd = "scp -r hkli@202.120.32.231:/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/data/%s/%s " \
#                   "%s/result/data/%s/"%(src_nm, src_sex, dd, src_nm,src_sex)
#             a = Popen(cmd, shell=True)
#             a.wait()
#             print(cmd)


datas = ["data_%d.hdf5"%i for i in range(10)]
for src_nm in sources:
    for dd in datas:
        cmd = "scp -r hkli@202.120.32.231:/mnt/ddnfs/data_users/hkli/selection_bias/paper_data/%s/result/data/%s " \
              "%s/result/data"%(src_nm, dd, src_nm)
        a = Popen(cmd, shell=True)
        a.wait()
        print(cmd)