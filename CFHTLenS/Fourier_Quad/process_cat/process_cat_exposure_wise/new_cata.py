from sys import path
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py
import time
import os


with open("/home/hklee/work/CFHT/CFHTLens.tsv","r") as f:
    cc = f.readlines()

headers = cc[0].split(",")
print(headers)
head_list = ["#"]
for i, ii in enumerate(headers):
    if ii != "PZ_full":
        head_list.append(ii)
print(len(head_list))
headers = "\t".join(head_list)

pzs = []
paras = [headers]
print(paras)

for i in range(1, len(cc)):

    ic = cc[i].split("\n")[0].split(",\"")

    all_para_ = ic[0].split(",")
    field_id = all_para_[0]
    para_ = all_para_[1:]

    pz_, mag = ic[1].split("\"")
    pz_ = pz_.split(",")
    pzs.append("\t".join(pz_) + "\n")

    all_para_[0] = all_para_[0].replace("W", "1")
    if "p" in all_para_[0]:
        all_para_[0] = all_para_[0].replace("p", "2")
    if "m" in all_para_[0]:
        all_para_[0] = all_para_[0].replace("m", "3")
    all_para_[0] = all_para_[0].replace("_", "")

    sub_paras = []
    sub_paras.extend(all_para_)
    sub_paras.append(mag.split(",")[1])
    paras.append("\t".join(sub_paras) + "\n")
    # print(sub_paras, len(sub_paras), "\n\n")

with open("/home/hklee/work/CFHT/CFHTLens.cat","w") as f:
    f.writelines(paras)
with open("/home/hklee/work/CFHT/CFHTLens_pz.cat","w") as f:
    f.writelines(pzs)