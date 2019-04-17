import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
from Fourier_Quad import Fourier_Quad
import tool_box
import h5py
from mpi4py import MPI
import numpy
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

area_num = 4
area_id = 1

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

block_scale = 60
margin = 0.5

flux_alt_thresh = 8.77
nstar_thresh = 12
total_area_thresh = 7
field_g1_bound = 0.005
field_g2_bound = 0.0075
c2_correction = 0.000498

fq = Fourier_Quad(64, 123)


h5f = h5py.File(data_path+"cata_result_ext.hdf5","r")
cat_data = h5f["/w_%d"%area_id].value
h5f.close()

# cut off
flux_alt_idx = cat_data[:, flux_alt_lb] >= flux_alt_thresh
nstar_idx = cat_data[:, nstar_lb] >= nstar_thresh
total_area_idx = cat_data[:, total_area_lb] >= total_area_thresh

fg1 = numpy.abs(cat_data[:, field_g1_lb])
fg2 = numpy.abs(cat_data[:, field_g2_lb])

fg1_idx = fg1 <= field_g1_bound
fg2_idx = fg2 <= field_g2_bound

cut_idx = flux_alt_idx & nstar_idx & total_area_idx & fg1_idx & fg2_idx

fg1 = fg1[cut_idx]
fg2 = fg2[cut_idx]

# the esitmators
mn = cat_data[:, mn_lb][cut_idx]
mu = cat_data[:, mu_lb][cut_idx]
mv = cat_data[:, mv_lb][cut_idx]
mg1 = cat_data[:, mg1_lb][cut_idx] - fg1*(mn+mu) - fg2*mv
mg2 = cat_data[:, mg2_lb][cut_idx] - (fg2 + c2_correction)*(mn-mu) - fg1*mv

ra = cat_data[:,ra_lb][cut_idx] * 60
dec = cat_data[:,dec_lb][cut_idx] * 60

ra_min, ra_max = ra.min() - margin, ra.max() + margin
dec_min, dec_max = dec.min() - margin, dec.max() + margin

grid_rows = int((dec_max - dec_min) / block_scale + 1)
grid_cols = int((ra_max - ra_min) / block_scale + 1)

dec_bin = numpy.array([dec_min + i * block_scale for i in range(grid_rows + 1)])
ra_bin = numpy.array([ra_min + i * block_scale for i in range(grid_cols + 1)])

if rank == 0:

    fig, ax = plt.subplots(figsize=(10, 10))
    nums = numpy.histogram2d(dec, ra, [dec_bin, ra_bin])[0]
    numpy.savez(pic_path+"w_%d_num.npz"%(rank+1), nums)
    im = ax.imshow(nums)
    fig.colorbar(im)
    plt.savefig(pic_path+"w_%d_num.png"%(rank+1),bbox_inches='tight')
    plt.close()
# [g1, g2]
shear_array = numpy.zeros((grid_rows, grid_cols*4))
grid_ids = [(i,j) for i in range(grid_rows) for j in range(grid_cols)]
my_grid = tool_box.allot(grid_ids, cpus)[rank]

for grid in my_grid:
    row, col = grid
    dec_1, dec_2 = dec_min + row*block_scale, dec_min + (row + 1)*block_scale
    ra_1, ra_2 = ra_min + col * block_scale, ra_min + (col + 1) * block_scale

    idx_r1 = dec >= dec_1
    idx_r2 = dec < dec_2

    idx_c1 = ra >= ra_1
    idx_c2 = ra < ra_2

    idxs = idx_r1 & idx_r2 & idx_c1 & idx_c2

    mg1_ = mg1[idxs]
    mg2_ = mg2[idxs]
    mnu1 = mn[idxs] + mu[idxs]
    mnu2 = mn[idxs] - mu[idxs]
    print(row, col, len(mg1_))
    g1, g1_sig = fq.fmin_g_new(mg1_, mnu1, bin_num=8)
    g2, g2_sig = fq.fmin_g_new(mg2_, mnu2, bin_num=8)

    shear_array[row, col] = g1
    shear_array[row, col+grid_cols] = g1_sig
    shear_array[row, col+grid_cols*2] = g2
    shear_array[row, col+grid_cols*3] = g2_sig

if rank > 0:
    comm.Send(shear_array, dest=0, tag=rank)
else:
    shear_final = numpy.zeros_like(shear_array)
    shear_final += shear_array
    for procs in range(1, cpus):
        recvs = numpy.empty((grid_rows, grid_cols*4), dtype=numpy.double)
        comm.Recv(recvs, source=procs, tag=procs)
        shear_final += recvs
        shear_array = numpy.row_stack((shear_array, recvs))
    numpy.savez(data_path + "w_%d_shear.npz"%area_id,
                shear_final, shear_array,numpy.array([grid_rows, grid_cols]), numpy.array([ra_min, ra_max, dec_min, dec_max]))

