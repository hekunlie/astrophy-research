import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
import matplotlib
import matplotlib.pyplot as plt
import plot_tool
import h5py
import numpy
from mpi4py import MPI

#
# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# cpus = comm.Get_size()

rank = 0

area_id = 1

# data_path = "/mnt/perc/hklee/CFHT/gg_lensing/data/"
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/"

h5f_pre = h5py.File(data_path+"cata_result_ext_cut.hdf5","r")
h5f_grid = h5py.File(data_path+ "cata_result_ext_grid.hdf5","r")


data_name = ["Z", "DISTANCE", "RA", "DEC", "COS_DEC", "G1", "G2", "N", "U", "V",
             "num_in_block", "block_start", "block_end", "block_boundy", "block_boundx"]

z_pre = h5f_pre["/w_%d/Z"%area_id].value
dist_pre = h5f_pre["/w_%d/DISTANCE"%area_id].value

ra_pre = h5f_pre["/w_%d/RA"%area_id].value
dec_pre = h5f_pre["/w_%d/DEC"%area_id].value
cos_dec_pre = h5f_pre["/w_%d/COS_DEC"%area_id].value

g1_pre = h5f_pre["/w_%d/G1"%area_id].value
g2_pre = h5f_pre["/w_%d/G2"%area_id].value
n_pre = h5f_pre["/w_%d/N"%area_id].value
u_pre = h5f_pre["/w_%d/U"%area_id].value
v_pre = h5f_pre["/w_%d/V"%area_id].value


z_back = h5f_grid["/background/w_%d/Z"%area_id].value
dist_back = h5f_grid["/background/w_%d/DISTANCE"%area_id].value

ra_back = h5f_grid["/background/w_%d/RA"%area_id].value
dec_back = h5f_grid["/background/w_%d/DEC"%area_id].value
cos_dec_back = h5f_grid["/background/w_%d/COS_DEC"%area_id].value

g1_back = h5f_grid["/background/w_%d/G1"%area_id].value
g2_back = h5f_grid["/background/w_%d/G2"%area_id].value
n_back = h5f_grid["/background/w_%d/N"%area_id].value
u_back = h5f_grid["/background/w_%d/U"%area_id].value
v_back = h5f_grid["/background/w_%d/V"%area_id].value

num_in_block = h5f_grid["/background/w_%d/num_in_block"%area_id].value
block_start = h5f_grid["/background/w_%d/block_start"%area_id].value
block_end = h5f_grid["/background/w_%d/block_end"%area_id].value
block_boundy = h5f_grid["/background/w_%d/block_boundy"%area_id].value
block_boundx = h5f_grid["/background/w_%d/block_boundx"%area_id].value

z_fore = h5f_grid["/foreground/w_%d/Z"%area_id].value
dist_fore = h5f_grid["/foreground/w_%d/DISTANCE"%area_id].value

ra_fore = h5f_grid["/foreground/w_%d/RA"%area_id].value
dec_fore = h5f_grid["/foreground/w_%d/DEC"%area_id].value
cos_dec_fore = h5f_grid["/foreground/w_%d/COS_DEC"%area_id].value

h5f_pre.close()
h5f_grid.close()

ra_min, ra_max = ra_back.min(), ra_back.max()
dec_min, dec_max = dec_back.min(), dec_back.max()
print("Background:")
print("RA: %6.5f ~ %6.5f"%(ra_min, ra_max))
print("DEC: %6.5f ~ %6.5f\n"%(dec_min, dec_max))

speed_c = 2.99792458 * 1e5
H_0 = 100
coeff = numpy.pi / 180 * speed_c / H_0
count = 0

G_t_crits = []
G_x_crits = []
crits = []
foregal_num = 300#z_fore.shape[0]

# m, n = divmod(foregal_num, cpus)

my_gal_s = 0#rank*m
my_gal_e = foregal_num#(rank+1)*m
# if rank < cpus - 1:
#     my_gal_e = my_gal_e + n

radius = [0.04, 0.06]

for i in range(my_gal_s, my_gal_e):
    npt = i
    ra_p, dec_p, cos_dec_p, z_p, dist_p = ra_fore[npt], dec_fore[npt], cos_dec_fore[npt], z_fore[npt], dist_fore[npt]
    dist_len = dist_p * speed_c / H_0
    #     print("The foreground galaxy.")
    #     print("RA: %7.5f, DEC: %7.5f, Cos_DEC: %7.5f"%(ra_p, dec_p, cos_dec_p))
    #     print("Redshift: %8.7f."%z_p)
    #     print("Distance: %8.7f*c/H_0. = %8.7f Mpc/h"%(dist_p, dist_len))
    #     print(coeff)

    idx1 = ra_back >= ra_p - 1
    idx2 = ra_back <= ra_p + 1

    idx3 = dec_back >= dec_p - 1
    idx4 = dec_back <= dec_p + 1

    idx_s = idx1 & idx2 & idx3 & idx4
    #     img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
    #     img.plot_img(1,1)
    #     img.axs[0][0].scatter(ra_back[idx_s], dec_back[idx_s])

    delta_ra = (ra_back[idx_s] - ra_p) * cos_dec_back[idx_s]
    delta_dec = dec_back[idx_s] - dec_p
    delta_theta = numpy.sqrt(delta_ra ** 2 + delta_dec ** 2)

    # idx_theta = delta_theta < 0.5
    # img.axs[0][0].scatter(ra_back[idx_s][idx_theta], dec_back[idx_s][idx_theta])

    diff_r = delta_theta * dist_back[idx_s] * coeff

    diff_z = z_back[idx_s] - z_p

    idx_z = diff_z > 0.1

    idx_r1 = diff_r >= radius[0]
    idx_r2 = diff_r < radius[1]

    idx_r = idx_r1 & idx_r2 & idx_z

    if idx_r.sum() > 0:
        #          rotation

        cos_2theta = (delta_ra[idx_r] / delta_theta[idx_r]) ** 2 - (delta_dec[idx_r] / delta_theta[idx_r]) ** 2
        sin_2theta = (delta_ra[idx_r] / delta_theta[idx_r]) * (delta_dec[idx_r] / delta_theta[idx_r]) * 2

        G_t = g1_back * cos_2theta - g2_back * sin_2theta

        G_x = g1_back * sin_2theta + g2_back * cos_2theta

        #     the comoving distance of source, the integrate without c/H_0
        dist_s = dist_back[idx_s][idx_r]

        crit_ = dist_s / (dist_s - dist_p) / dist_p / (1 + z_p)

        G_t_crit = (G_t * crit_).tolist()
        G_x_crit = (G_x * crit_).tolist()

        G_t_crits.extend(G_t_crit)
        G_x_crits.extend(G_x_crit)
        crits.extend(crit_.tolist())
        count += idx_r.sum()
        # print("Crit: ",crit_, "G_t",G_t_crit,"G_x", G_x_crit)
        print("%d galaxies found in [%.4f, %.4f] Mpc \n" % (idx_r.sum(), radius[0], radius[1]))
print(G_t_crits)
print(G_x_crits)
print(crits)
print(count,len(G_t_crits), len(G_x_crits), len(crits))
numpy.savez("data.npz",G_t_crits,G_x_crits,crits)
# G_t_crits = comm.gather(G_t_crits, root=0)
# G_x_crits = comm.gather(G_x_crits, root=0)
# crits = comm.gather(crits, root=0)
#
# if rank == 0:
#     G_t_total = []
#     G_x_total = []
#     crits_total = []
#
#     for i in range(cpus):
#         G_t_total.extend(G_t_crits[i])
#         G_x_total.extend(G_x_crits[i])
#         crits_total.extend(crits[i])
#
#     print(G_t_total)
#     print(G_x_total)
#     print(crits_total)

# # print(diff_r)
# # print(diff_r[idx_r])
# img.axs[0][0].scatter(ra_back[idx_s][idx_r], dec_back[idx_s][idx_r])
# img.axs[0][0].scatter(ra_p, dec_p, s=15)

# dx = 0.2
# img.axs[0][0].set_ylim(dec_p- dx, dec_p+dx)
# img.axs[0][0].set_xlim(ra_p- dx, ra_p+dx)