import numpy

size = 48
cen = int((size * size + size) / 2)  # center of image

for i in range(1):
    for j in range(1):
        arr = numpy.zeros((size, size))
        arr[int(size / 2), int(size / 2)] = 0.5
        pos = []
        tag = 0
        pk = 0
        z = []
        x_d = []
        my, mx = numpy.mgrid[-2:3, -2:3]
        x, y = mx.reshape((1, 25)), my.reshape((1, 25))
        for ix in [0, 4, 20, 24]:
            x[0,ix] = 0
            y[0,ix] = 0
        for m in range(-2, 3):
            p = (i + m + size) % size
            for n in range(-2, 3):
                q = (j + n + size) % size
                if p * size + q == cen:
                    pk = tag
                    x[0,tag] = 0
                    y[0,tag] = 0
                tag += 1
        if pk != 0:
            num =