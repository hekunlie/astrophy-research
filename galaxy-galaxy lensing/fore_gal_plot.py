import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import tool_box
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

pic_path = data_path + "pic/"

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

block_scale = 40
margin = 0.5

# ra dec flag  flux_radius e1 e2 weight  fitclass  SNR  MASK  Z_B  m  c2 LP_Mi  star_flag  MAG_i
# 0   1   2     3          4  5  6       7         8    9     10   11 12 13      14        15
h5f = h5py.File(cfht_cata_path+"cata.hdf5","r")
data = h5f["/w_%d"%(rank+1)].value
h5f.close()

ra = data[:, 0]*60
dec = data[:, 1]*60

mask = data[:, 9]
mask.astype(dtype=numpy.intc)

redshift = data[:, 10]
idx1 = redshift > 0

star_flag = data[:, 14]
star_flag.astype(dtype=numpy.intc)

mag = data[:, 15]
idx2 = mag > 0

# print(rank, mask.shape[0],mag.min(), mag.max(), redshift.min(), redshift.max(),idx1.sum(), idx2.sum())
# the grid: x-axis--ra, y-axis--dec
ra_min, ra_max = ra.min() - margin, ra.max() + margin
dec_min, dec_max = dec.min() - margin, dec.max() + margin

# the grid
grid_rows = int((dec_max - dec_min) / block_scale + 1)
grid_cols = int((ra_max - ra_min) / block_scale + 1)

# the foreground
fig, ax = plt.subplots(figsize=(15,15))
idx_m1 = mag < 19.5
idx_m2 = mag > 0

idx_z = redshift < 0.3

idxs = idx_m1&idx_m2&idx_z

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

plt.savefig(pic_path+"w_%d_fore.png"%(rank+1))
plt.close()
print("fore",len(ra[idxs]))

# all galaxies
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

idx_m1 = mag < 30
idx_m2 = mag > 0
idxs = idx_m1&idx_m2
print(rank, mag[idxs].max(), mag[idxs].min())
ax.scatter(x=ra[idxs], y=dec[idxs],s=5, cmap=cmap)
plt.savefig(pic_path+"w_%d_all.png"%(rank+1))
plt.close()

# number count in ra-dec grid
fig, ax = plt.subplots(figsize=(10, 10))
idx_m1 = mag < 25
idx_m2 = mag > 0
idx_g = star_flag == 0
idxs = idx_m1&idx_m2&idx_g

# print(rank, mag.shape[0], idxs.sum())
h5f = h5py.File(data_path+"cata_result_ext.hdf5", "r")
data = h5f["/w_%d"%(rank+1)].value
h5f.close()

ra = data[:,ra_lb]*60
dec = data[:,dec_lb]*60
dec_bin = numpy.array([dec_min + i * block_scale for i in range(grid_rows + 1)])
ra_bin = numpy.array([ra_min + i * block_scale for i in range(grid_cols + 1)])
# y,x,y_bin,x_bin
nums = numpy.histogram2d(dec, ra, [dec_bin, ra_bin])[0]
numpy.savez(pic_path+"w_%d_num.npz"%(rank+1), nums)
inv = range(nums.shape[0]-1, -1, -1)
im = ax.imshow(nums[inv])
fig.colorbar(im)
plt.savefig(pic_path+"w_%d_all_bin.png"%(rank+1),bbox_inches='tight')
plt.close()