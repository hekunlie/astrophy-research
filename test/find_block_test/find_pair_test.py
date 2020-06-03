import numpy
import matplotlib.pyplot as plt
import time

seed = 112312

rng = numpy.random.RandomState(seed)
# source
src_num = 10000
src_x = rng.uniform(0, 2400, 10000)
src_y = rng.uniform(0, 4000, 10000)

# set up bins for source xy
grid_dx, grid_dy = 5,5
nx = int((src_x.max()-src_x.min())/grid_dx)+1
ny = int((src_y.max()-src_y.min())/grid_dy)+1
x_min, x_max = src_x.min(), src_x.max()
y_min, y_max = src_y.min(), src_y.max()
x_bin = numpy.array([x_min + grid_dx*i for i in range(nx+1)])
y_bin = numpy.array([y_min + grid_dy*i for i in range(ny+1)])

grid_label = numpy.zeros((src_num,),dtype=numpy.intc)-1

grid_cent = numpy.zeros((ny*nx, 2))

x_grid = numpy.zeros((ny*nx, 4))
y_grid = numpy.zeros((ny*nx, 4))

for i in range(ny):
    idx1 = src_y >= y_bin[i]
    idx2 = src_y < y_bin[i+1]
    idx_ = idx1&idx2
    for j in range(nx):
        idx3 = src_x >= x_bin[j]
        idx4 = src_x < x_bin[j+1]
        idx = idx3&idx4&idx_

        tag = i*nx + j

        grid_label[idx] = tag

        grid_cent[tag,0] = (y_bin[i] + y_bin[i+1])/2
        grid_cent[tag,1] = (x_bin[j] + x_bin[j+1])/2

        x_grid[tag,0] = x_min + grid_dx*i
        x_grid[tag,1] = x_min + grid_dx*(i+1)
        x_grid[tag,2] = x_min + grid_dx*(i+1)
        x_grid[tag,3] = x_min + grid_dx*i

        y_grid[tag,0] = y_min + grid_dy*j
        y_grid[tag,1] = y_min + grid_dy*j
        y_grid[tag,2] = y_min + grid_dy*(j+1)
        y_grid[tag,3] = y_min + grid_dy*(j+1)

print(grid_label[grid_label<0])

fore_num = 100
fore_x = rng.uniform(-200, -100, fore_num)
fore_y = rng.uniform(-200, -100, fore_num)
pair_num = numpy.zeros((fore_num, ))

plt.scatter(src_x, src_y)
plt.scatter(fore_x, fore_y)
plt.show()

rad = [300, 1000]
rad_sq = [300**2, 1000**2]
loops = 5

t1 = time.time()
for k in range(loops):
    for i in range(fore_num):
        # simulate the process in C++
        for j in range(src_num):
            radius = numpy.sqrt((src_x[j] - fore_x[i])**2 + (src_y[j] - fore_y[i])**2)
            if radius >= rad[0] and radius < rad[1]:
                pair_num[i] += 1


        # radius = (src_x - fore_x[i]) ** 2 + (src_y - fore_y[i]) ** 2
        # idx1 = radius >= rad[0]
        # idx2 = radius < rad[1]
        # idx = idx1&idx2
        # pair_num[i] = idx.sum()
t2 = time.time()
print(t2-t1)
print(pair_num)

t3 = time.time()
for k in range(loops):
    for i in range(fore_num):
        # simulate the process in C++
        for j in range(src_num):
            radius = (src_x[j] - fore_x[i])**2 + (src_y[j] - fore_y[i])**2
            if radius >= rad_sq[0] and radius < rad_sq[1]:
                pair_num[i] += 1

        # radius = (src_x - fore_x[i]) ** 2 + (src_y - fore_y[i]) ** 2
        # idx1 = radius >= rad[0]
        # idx2 = radius < rad[1]
        # idx = idx1&idx2
        # pair_num[i] = idx.sum()
t4 = time.time()
print(t4-t3)
print(pair_num)



exit()
t3 = time.time()
for k in range(loops):
    for i in range(fore_num):
        # simulate the process in C++
        for j in range(src_num):
            radius = (src_x[j] - fore_x[i])**2 + (src_y[j] - fore_y[i])**2
            if radius >= rad[0] and radius < rad[1]:
                pair_num[i] += 1

t4 = time.time()
print(t3-t4)
print(pair_num)