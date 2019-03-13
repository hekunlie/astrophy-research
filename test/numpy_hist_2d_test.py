import numpy
import matplotlib.pyplot as plt


num = 100
block_scale = 4
dec_min, dec_max = -5, 20
ra_min, ra_max = -10, 20
grid_rows = int((dec_max - dec_min) / block_scale + 1)
grid_cols = int((ra_max - ra_min) / block_scale + 1)
ra = numpy.random.random_sample(num)*(ra_max - ra_min)+ra_min
dec = numpy.random.random_sample(num)*(dec_max - dec_min)+dec_min

dec_bin = numpy.array([dec_min + i * block_scale for i in range(grid_rows + 1)])
ra_bin = numpy.array([ra_min + i * block_scale for i in range(grid_cols + 1)])
# y, x ,[y_bin, x_bin]
nums = numpy.histogram2d(dec, ra, [dec_bin, ra_bin])[0]
print(numpy.sort(ra))
print(ra_bin)
print(numpy.sort(dec))
print(dec_bin)

fig = plt.figure(figsize=(12,6))
ax1 = fig.add_subplot(121)
ax1.imshow(nums)

ax2 = fig.add_subplot(122)
for i in range(grid_rows + 1):
    # ax.plot([x1, x2..],[y1, y2,..])
    # the horizontal line
    ax2.plot([ra_min, ra_min + grid_cols * block_scale],
            [dec_min + i * block_scale, dec_min + i * block_scale], c="black", linestyle="--", linewidth=0.8)
for j in range(grid_cols + 1):
    # the vertical line
    ax2.plot([ra_min + j * block_scale, ra_min + j * block_scale],
            [dec_min, dec_min + grid_rows * block_scale], c="black", linestyle="--", linewidth=0.8)
ax2.scatter(ra,dec,s=6)
plt.show()