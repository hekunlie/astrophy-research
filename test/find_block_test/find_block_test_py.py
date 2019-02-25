from sys import path
path.append('D:/GitHub/astrophy-research/my_lib')
import time
import tool_box
import numpy
import matplotlib.pyplot as plt

scale = 5
ny,nx = 10,10
radius_s, radius_e = 7.5, 38.7

# grid
y_zero, x_zero = -100, 50
my,mx = numpy.mgrid[0:ny+1,0:nx+1]*scale
my += y_zero
mx += x_zero

# boundary
boundx = numpy.zeros((ny*nx, 4))
boundy = numpy.zeros((ny*nx, 4))
for i in range(ny):
    for j in range(nx):
        squen = i*nx + j
        boundy[squen] = my[i, j], my[i, j+1], my[i+1, j], my[i+1, j+1]
        boundx[squen] = mx[i, j], mx[i, j+1], mx[i+1, j], mx[i+1, j+1]


theta = numpy.linspace(0,2*numpy.pi,5000)
for iy in range(ny):
    for ix in range(nx):
        plt.figure(figsize=(8, 8))
        plt.axes().set_aspect('equal', 'datalim')
        for i in range(ny + 1):

            plt.plot([x_zero, x_zero+nx * scale], [y_zero+i * scale, y_zero+i * scale], c="black", linewidth=0.5)
            for j in range(nx + 1):
                plt.plot([x_zero+j * scale, x_zero+j * scale], [y_zero, y_zero+ny * scale], c="black", linewidth=0.5)
        dec = y_zero+iy * scale + numpy.random.rand() * scale
        ra = x_zero+ix * scale + numpy.random.rand() * scale
        # dec, ra = 104.68487, 54.31511
        plt.scatter(ra,dec,s=1)
        plt.title("%.4f,%.4f"%(dec - iy*scale, ra-ix*scale))
        target_blocks = tool_box.find_block(scale, radius_s, radius_e, iy, ix, dec, ra,
                                                    ny, nx, boundy, boundx)
        # show the found blocks
        for blks in target_blocks:
            blk_x, blks_y = x_zero+(blks[1]+0.5)*scale,y_zero+(blks[0]+0.5)*scale
            plt.scatter(blk_x, blks_y, s=5,c="blue")

        rs_y, rs_x = dec + radius_s*numpy.sin(theta), ra+radius_s*numpy.cos(theta)
        sy1 = rs_y >= y_zero
        sy2 = rs_y <= y_zero+ny*scale
        sx1 = rs_x >= x_zero
        sx2 = rs_x <= x_zero+nx*scale
        idx1 = sy1&sy2&sx1&sx2
        re_y, re_x = dec + radius_e * numpy.sin(theta), ra + radius_e * numpy.cos(theta)
        ey1 = re_y >= y_zero
        ey2 = re_y <= y_zero + ny*scale
        ex1 = re_x >= x_zero
        ex2 = re_x <= x_zero+nx*scale
        idx2 = ey1&ey2&ex1&ex2
        plt.scatter(rs_x[idx1], rs_y[idx1], c="coral",s=0.1)
        plt.scatter(re_x[idx2], re_y[idx2], c="coral", s=0.1)
        print(target_blocks)
        plt.show()
        plt.close()