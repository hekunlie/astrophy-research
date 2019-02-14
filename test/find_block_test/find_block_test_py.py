from sys import path
path.append('E:/GitHub/astrophy-research/my_lib')
import time
import tool_box
import numpy
import matplotlib.pyplot as plt

scale = 5
ny,nx = 10,10
radius_s, radius_e = 7.5, 38.7

# grid
my,mx = numpy.mgrid[0:ny+1,0:nx+1]*scale

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

            plt.plot([0, nx * scale], [i * scale, i * scale], c="black", linewidth=0.5)
            for j in range(nx + 1):
                plt.plot([j * scale, j * scale], [0, ny * scale], c="black", linewidth=0.5)
        # dec = iy * scale + numpy.random.rand() * scale
        # ra = ix * scale + numpy.random.rand() * scale
        dec, ra = 4.68487, 4.31511
        plt.scatter(ra,dec,s=1)
        plt.title("%.4f,%.4f"%(dec - iy*scale, ra-ix*scale))
        target_blocks = tool_box.find_block(scale, radius_s, radius_e, iy, ix, dec, ra,
                                                    ny, nx, boundy, boundx)
        # show the found blocks
        for blks in target_blocks:
            blk_x, blks_y = (blks[1]+0.5)*scale,(blks[0]+0.5)*scale
            plt.scatter(blk_x, blks_y, s=5,c="blue")

        rs_y, rs_x = dec + radius_s*numpy.sin(theta), ra+radius_s*numpy.cos(theta)
        sy1 = rs_y >= 0
        sy2 = rs_y <= ny*scale
        sx1 = rs_x >=0
        sx2 = rs_x <= nx*scale
        idx1 = sy1&sy2&sx1&sx2
        re_y, re_x = dec + radius_e * numpy.sin(theta), ra + radius_e * numpy.cos(theta)
        ey1 = re_y >= 0
        ey2 = re_y <= ny*scale
        ex1 = re_x >=0
        ex2 = re_x <= nx*scale
        idx2 = ey1&ey2&ex1&ex2
        plt.scatter(rs_x[idx1], rs_y[idx1], c="coral",s=0.1)
        plt.scatter(re_x[idx2], re_y[idx2], c="coral", s=0.1)
        print(target_blocks)
        plt.show()
        plt.close()