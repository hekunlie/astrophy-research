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
print(headers)

pzs = []
paras = [headers]
for i in range(1,len(cc)):
#     print(cc[i])
    ic = cc[i].split("\n")[0].split(",\"")
    para_ = ic[0].split(",")[1:]
#     print(para_)
    pz_, mag = ic[1].split("\"")
    pz_ = pz_.split(",")
    para_.append(mag.split(",")[1])
#     print(pz_, mag)
#     print(len(pz_))
#     print(para_,len(para_))
    paras.append("\t".join(para_) + "\n")
    pzs.append("\t".join(pz_)+ "\n")
# print(paras)
# print(pzs)

with open("/home/hklee/work/CFHT/CFHTLens.cat","w") as f:
    f.writelines(paras)
with open("/home/hklee/work/CFHT/CFHTLens_pz.cat","w") as f:
    f.writelines(pzs)