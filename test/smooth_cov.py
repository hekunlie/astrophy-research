import numpy
import matplotlib.pyplot as plt
import h5py

size = 48
cen = int((size * size + size) / 2)  # center of image

tag = 0
all_matrix = numpy.zeros((25*6, 25))
for m in range(0, 5): # y
    for n in range(0, 5): # x
        y, x = numpy.mgrid[-2:3, -2:3]
        for iy in [0, 4]:
            for ix in [0, 4]:
                x[iy, ix] = 0
                y[iy, ix] = 0
        if tag in [0, 4, 20, 24]:
            num = 21
        else:
            num = 20
        x[m, n] = 0
        y[m, n] = 0

        x1 = numpy.sum(x)
        x2 = numpy.sum(x*x)
        x3 = numpy.sum(x*x*x)
        x4 = numpy.sum(x*x*x*x)
        y1 = numpy.sum(y)
        y2 = numpy.sum(y*y)
        y3 = numpy.sum(y*y*y)
        y4 = numpy.sum(y*y*y*y)
        xy = numpy.sum(x*y)
        x2y = numpy.sum(x*x*y)
        xy2 = numpy.sum(x*y*y)
        x3y = numpy.sum(x*x*x*y)
        x2y2 = numpy.sum(x*x*y*y)
        xy3 = numpy.sum(x*y*y*y)

        cov = numpy.array([[num, x1, y1, x2, xy, y2],
                           [x1, x2, xy, x3, x2y, xy2],
                           [y1, xy, y2, x2y, xy2, y3],
                           [x2, x3, x2y, x4, x3y, x2y2],
                           [xy, x2y, xy2, x3y, x2y2, xy3],
                           [y2, xy2, y3, x2y2, xy3, y4]])
        inv_cov = numpy.linalg.inv(cov)
        print(inv_cov[0])
        sub_matrix = numpy.ones((6,25))*inv_cov[0,0] # the first row are const.
        sub_matrix[1] = inv_cov[0, 1] * (x.reshape((1, 25)))
        sub_matrix[2] = inv_cov[0, 2] * (y.reshape((1, 25)))
        sub_matrix[3] = inv_cov[0, 3] * ((x*x).reshape((1, 25)))
        sub_matrix[4] = inv_cov[0, 4] * ((x*y).reshape((1, 25)))
        sub_matrix[5] = inv_cov[0, 5] * ((y*y).reshape((1, 25)))
        all_matrix[tag*6: (tag+1)*6] = sub_matrix
        print(sub_matrix)
        # print(all_matrix)
        # print(x)
        # print(y)
        # print(cov)
        #
        # plt.subplot(121)
        # plt.imshow(x)
        # plt.subplot(122)
        # plt.title(str(num))
        # plt.imshow(y)
        # plt.show()

        tag += 1
f = h5py.File("coeffs.hdf5","w")
f["/data"] = all_matrix
f.close()