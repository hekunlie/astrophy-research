import pandas as pd
import numpy as np
import os
import time
t1=time.clock()
cor_data = np.zeros((10000,2))
df1 = pd.read_excel("E:/0_chip_13.xlsx")
data = df1.values
col = ['BJ CORRECTION','RE CORRECTION']
cor_data[:, 0] = 2-(data[:, 1]**2+data[:, 7]**2)
cor_data[:, 1] = 2-(data[:, 2]**2+data[:, 8]**2)
df = pd.DataFrame(data=cor_data,columns=col)
new_df = pd.concat([df1,df],axis=1)
new_df.to_excel("E:/0_chip_13.xlsx")

t2=time.clock()
print(t2-t1)