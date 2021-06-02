from sys import path
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py
import time
import os

zbin = numpy.array([0.025 + i * 0.05 for i in range(70)])

with open("/home/hklee/work/CFHT/CFHTLens.tsv","r") as f:
    cc = f.readlines()

src_num = len(cc) - 1

headers = cc[0].split("\n")[0].split(",")
print(headers)
head_list = ["#"]
for i, ii in enumerate(headers):
    if ii != "PZ_full":
        head_list.append(ii)
head_list.append("Z_E\n")

print(len(head_list))
headers = "\t".join(head_list)

pzs = []
paras = [headers]
print(paras)

pz_arr = numpy.zeros((src_num, 70),dtype=numpy.float32)

for i in range(1, len(cc)):

    ic = cc[i].split("\n")[0].split(",\"")

    all_para_ = ic[0].split(",")
    field_id = all_para_[0]
    para_ = all_para_[1:]

    pz_, mag = ic[1].split("\"")
    pz_ = pz_.split(",")
    for ii in range(len(pz_)):
        pz_arr[i-1, ii] = float(pz_[ii])

    pz_sum = numpy.sum(pz_arr[i-1])
    if pz_sum != 0:
        ze = numpy.sum(zbin*pz_arr[i-1])/pz_sum
    else:
        ze = -99.0

    all_para_[0] = all_para_[0].replace("W", "1")
    if "p" in all_para_[0]:
        all_para_[0] = all_para_[0].replace("p", "2")
    if "m" in all_para_[0]:
        all_para_[0] = all_para_[0].replace("m", "3")
    all_para_[0] = all_para_[0].replace("_", "")

    sub_paras = []
    sub_paras.extend(all_para_)
    sub_paras.append(mag.split(",")[1])
    sub_paras.append("%f"%ze)
    paras.append("\t".join(sub_paras) + "\n")
    # print(sub_paras, len(sub_paras), "\n\n")

with open("/home/hklee/work/CFHT/CFHTLens.cat","w") as f:
    f.writelines(paras)

numpy.savez("/home/hklee/work/CFHT/CFHTLens_pz.npz", pz_arr)



