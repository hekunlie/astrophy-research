import numpy
import time
import os
import pandas
files = os.listdir('/lmc/selection_bias/result/data/')
paths = []
path = '/lmc/selection_bias/result/data/'
gg1 = []
gg2 = []
for i in files:
    if ".xlsx" in i:
        df = pandas.read_excel(path+i).values
        g1 = df[1,5]
        g2 = df[1,11]
        if g1 not in gg1:
            gg1.append(g1)
        if g2 not in gg2:
            gg2.append(g2)

print(gg1)
print(gg2)