import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
import h5py
from mpi4py import MPI
import numpy
import h5py
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

area_num = 4

result_source = "result_ext"

envs_path = "%s/work/envs/envs.dat"%my_home

gets_item = [["cfht", "cfht_path_catalog", "0"], ["gg_lensing", "ggl_path_data", "0"]]
path_items = tool_box.config(envs_path, ["get", "get"], gets_item)

cata_path, data_path = path_items

cfht_cata_path = cata_path + "cfht_cata/"
fourier_cata_path = cata_path + "fourier_cata/"


gets_item = [["fresh_para_idx", "nstar", "0"], ["fresh_para_idx", "flux_alt", "0"],
             ["fresh_para_idx", "ra", "0"], ["fresh_para_idx", "dec", "0"],
             ["fresh_para_idx", "gf1", "0"], ["fresh_para_idx", "gf2", "0"],
             ["fresh_para_idx", "g1", "0"], ["fresh_para_idx", "g2", "0"],
             ["fresh_para_idx", "de", "0"], ["fresh_para_idx", "h1", "0"],
             ["fresh_para_idx", "h2", "0"], ["fresh_para_idx", "total_area", "0"]]
gets = ["get" for i in range(len(gets_item))]
para_items = tool_box.config(envs_path, gets, gets_item)

nstar_lb = int(para_items[0])
flux_alt_lb = int(para_items[1])
total_area_lb = int(para_items[11])

ra_lb = int(para_items[2])
dec_lb = int(para_items[3])

field_g1_lb = int(para_items[4])
field_g2_lb = int(para_items[5])

mg1_lb = int(para_items[6])
mg2_lb = int(para_items[7])
mn_lb = int(para_items[8])
mu_lb = int(para_items[9])
mv_lb = int(para_items[10])
z_lb = -3
mag_lb = -2

block_scale = 10
margin = 0.5
h5f = h5py.File(data_path+"cata_result_ext.hdf5","r")
data = h5f["/w_%d"%(rank+1)].value
h5f.close()

mag = data[:, mag_lb]
redshift = data[:, z_lb]
ra = data[:, ra_lb]*60
dec = data[:, dec_lb]*60

# the grid: x-axis--ra, y-axis--dec
ra_min, ra_max = ra.min() - margin, ra.max() + margin
dec_min, dec_max = dec.min() - margin, dec.max() + margin

# the grid
grid_rows = int((dec_max - dec_min) / block_scale + 1)
grid_cols = int((ra_max - ra_min) / block_scale + 1)


idx_m1 = mag < 20
idx_m2 = mag > 18.5

idx_z = redshift < 0.3

idxs = idx_m1&idx_m2&idx_z

fig, ax = plt.subplots(figsize=(15,15))

for i in range(grid_rows + 1):
    # ax.plot([x1, x2..],[y1, y2,..])
    # the horizontal line
    ax.plot([ra_min, ra_min + grid_cols * block_scale],
            [dec_min + i * block_scale, dec_min + i * block_scale], c="black", linestyle="--", linewidth=0.8)
for j in range(grid_cols + 1):
    # the vertical line
    ax.plot([ra_min + j * block_scale, ra_min + j * block_scale],
            [dec_min, dec_min + grid_rows * block_scale], c="black", linestyle="--", linewidth=0.8)

cmap = plt.cm.get_cmap('RdYlBu')
normalize = matplotlib.colors.Normalize(vmin=numpy.min(mag[idxs]), vmax=numpy.max(mag[idxs]))
colors = cmap(normalize(mag[idxs]))

ax.scatter(x=ra[idxs], y=dec[idxs],s=5, cmap=cmap)
cax, _ = matplotlib.colorbar.make_axes(ax)
cbar = matplotlib.colorbar.ColorbarBase(cax, cmap=cmap, norm=normalize)

plt.savefig(data_path+"w_%d_fore.png"%(rank+1))
plt.close()



fig, ax = plt.subplots(figsize=(15,15))

for i in range(grid_rows + 1):
    # ax.plot([x1, x2..],[y1, y2,..])
    # the horizontal line
    ax.plot([ra_min, ra_min + grid_cols * block_scale],
            [dec_min + i * block_scale, dec_min + i * block_scale], c="black", linestyle="--", linewidth=0.8)
for j in range(grid_cols + 1):
    # the vertical line
    ax.plot([ra_min + j * block_scale, ra_min + j * block_scale],
            [dec_min, dec_min + grid_rows * block_scale], c="black", linestyle="--", linewidth=0.8)

idx_m1 = mag < 3000
idx_m2 = mag > -1000

idxs = idx_m1&idx_m2&idx_z

ax.scatter(x=ra[idxs], y=dec[idxs],s=5, cmap=cmap)
plt.savefig(data_path+"w_%d_all.png"%(rank+1))
plt.close()


fig, ax = plt.subplots(figsize=(10,10))
dec_bin = numpy.array([dec_min + i * block_scale for i in range(grid_rows + 1)])
ra_bin = numpy.array([ra_min + i * block_scale for i in range(grid_cols + 1)])
nums = numpy.histogram2d(ra, dec, [ra_bin, dec_bin])[0]

im = ax.imshow(nums)
fig.colorbar(im)
plt.savefig(data_path+"w_%d_all_bin.png"%(rank+1))
plt.close()