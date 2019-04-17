import numpy
from sys import path
path.append("E:/Github/astrophy-research/mylib")
import matplotlib.pyplot as plt

data = numpy.load("E:/pic/w_1_num.npz")["arr_0"]
cluster = numpy.where(data == data.max())

arrs = numpy.load("E:/pic/w_1_shear.npz")
shear = arrs["arr_0"]
shear_all = arrs["arr_1"]
rc = arrs["arr_2"]
grid_rows, grid_cols = rc[0], rc[1]
ra_min, ra_max, dec_min, dec_max = arrs["arr_3"]
print(ra_min, ra_max, dec_min, dec_max)
block_scale = 60
print(shear.shape, grid_rows, grid_cols, cluster)

g1 = shear[:, :grid_cols]
g1_sig = shear[:, grid_cols:grid_cols * 2]
g2 = shear[:, grid_cols * 2:grid_cols * 3]
g2_sig = shear[:, grid_cols * 3:grid_cols * 4]

g = numpy.sqrt(g1*g1 + g2*g2)
g_max = g.max()

plt.subplot(231)
plt.title("g1_sig")
plt.imshow(g1_sig)
plt.subplot(232)
plt.title("g2_sig")
plt.imshow(g2_sig)
plt.subplot(233)
plt.title("g")
plt.imshow(g)

plt.subplot(234)
plt.title("g1")
plt.imshow(g1)
plt.subplot(235)
plt.title("g2")
plt.imshow(g2)
plt.subplot(236)
plt.title("Galaxy number")
plt.imshow(data)
plt.show()
plt.close()

fig, ax = plt.subplots(figsize=(10,10))

for i in range(grid_rows + 1):
    # ax.plot([x1, x2..],[y1, y2,..])
    # the horizontal line
    ax.plot([ra_min, ra_min + grid_cols * block_scale],
            [dec_min + i * block_scale, dec_min + i * block_scale], c="black", linestyle="--", linewidth=0.8)
for j in range(grid_cols + 1):
    # the vertical line
    ax.plot([ra_min + j * block_scale, ra_min + j * block_scale],
            [dec_min, dec_min + grid_rows * block_scale], c="black", linestyle="--", linewidth=0.8)

benchmark = block_scale
for i in range(grid_rows):
    for j in range(grid_cols):

        cen_x = j * block_scale + ra_min + 0.5 * block_scale
        cen_y = i * block_scale + dec_min + 0.5 * block_scale
        dx = (cluster[0][0] - j)*block_scale
        dy = (cluster[1][0] - i)*block_scale
        distance = numpy.sqrt(dx*dx+dy*dy)
        cos_phi = dx / distance
        sin_phi = dy / distance
        cos_2phi = cos_phi*cos_phi - sin_phi*sin_phi
        sin_2phi = 2*cos_phi*sin_phi

        g_t = -g1[i,j]*cos_2phi - g2[i,j]*sin_2phi
        g_x = g1[i,j]*sin_2phi - g2[i,j]*cos_2phi
        cos_phi = numpy.sqrt((1 + g_t / g_max) / 2)
        sin_phi = numpy.sqrt((1 - g_t / g_max) / 2)

        line_len = g[i, j]/g_max*benchmark
        cos_phi_ = numpy.sqrt((1 + g1[i, j] / g_max) / 2)
        sin_phi_ = numpy.sqrt((1 - g1[i, j] / g_max) / 2)
        dx_ = line_len/2*cos_phi
        dy_ = line_len/2*sin_phi

        plt.plot([cen_x-dx_, cen_x+dx_], [cen_y-dy_, cen_y+dy_],c="black", linewidth=1.6)
        pts_size = 3
        if i == cluster[0][0] and j == cluster[1][0]:
            pts_size = 20
        ax.scatter(cen_x, cen_y, s=pts_size, color="black")


plt.show()
plt.close()