import matplotlib.pyplot as plt
import numpy
import pandas


arr_s = numpy.loadtxt('E:/test.cat')
x = arr_s[:,-2]
y = arr_s[:,-1]
arr = numpy.zeros((100,100))
print(len(x))
for i in range(len(x)):
    mx,dx = divmod(x[i],60)
    my,dy = divmod(y[i],60)
    arr[int(my),int(mx)]+=10
    if mx==11 and my==40:
        print(y[i],x[i])
    # if arr[int(my),int(mx)]==20:
    #     print(y[i],x[i])
    #     print(my,mx)
plt.figure()
plt.imshow(arr)
plt.colorbar()
plt.show()