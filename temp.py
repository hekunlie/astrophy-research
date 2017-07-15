import numpy as np

a = [x for x in range(5)]
print(a,type(a),np.where(a==np.max(a))[0])