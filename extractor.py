from astropy.io import fits
import matplotlib.pyplot as plt
import numpy
import copy
my,mx = numpy.mgrid[0:50,0:50]
arr = numpy.exp(-(mx-15)**2/8-(my-15)**2/8)+numpy.exp(-(mx-35)**2/8-(my-35)**2/8)
idx = arr <0.1
arr[idx]=0
idx = arr>0
arr[idx]=1
plt.imshow(arr)
plt.colorbar()
plt.show()

def detect(mask,ini_x,ini_y):
    signal = []
    if mask[ini_x,ini_y]!=0:
        mask[ini_x,ini_y] = 0
        signal.extend([ini_x,ini_y])
        for cor in [(-1,0),(1,0),(0,-1),(0,1)]:
            detect(mask,ini_x+cor[0],ini_y+cor[1])
            return signal
signal =[]

mask = copy.copy(arr)
for i in range(50):
    for j in range(50):
        if arr[i,j]!=0:
            sig = detect(arr,i,j)
            signal.append(sig)

print(len(signal))
