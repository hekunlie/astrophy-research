import numpy
from astropy.io import fits
import h5py
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
import tool_box
from plot_tool import Image_Plot
import matplotlib.pyplot as plt


data_n = fits.open("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/boss/galaxy_DR12v5_CMASS_North.fits.gz")[1].data
data_s = fits.open("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/boss/galaxy_DR12v5_CMASS_South.fits.gz")[1].data

ra_n, dec_n, z_n = data_n["RA"], data_n["DEC"], data_n["Z"]

ra_n_min, ra_n_max = ra_n.min(), ra_n.max()
dec_n_min, dec_n_max = dec_n.min(), dec_n.max()
z_n_min, z_n_max = z_n.min(), z_n.max()

ra_s, dec_s, z_s = data_s["RA"], data_s["DEC"], data_s["Z"]

ra_s_min, ra_s_max = ra_s.min(), ra_s.max()
dec_s_min, dec_s_max = dec_s.min(), dec_s.max()
z_s_min, z_s_max = z_s.min(), z_s.max()


print("North: RA: %.3f ~ %.3f, DEC: %.3f ~ %.3f, Z: %.4f ~ %.4f"%(ra_n_min, ra_n_max,dec_n_min, dec_n_max,z_n_min, z_n_max))

print("South: RA: %.3f ~ %.3f, DEC: %.3f ~ %.3f, Z: %.4f ~ %.4f"%(ra_s_min, ra_s_max,dec_s_min, dec_s_max,z_s_min, z_s_max))

#  The original catalog
img = Image_Plot(fig_x=12, fig_y=9)
img.set_style()
img.plot_img(1,1)

img.axs[0][0].scatter(ra_n, dec_n, s=3, label="CMASS_NORTH")
img.axs[0][0].scatter(ra_s, dec_s, s=3, label="CMASS_SOUTH")

h5f = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/catalog/cfht_cata/cata.hdf5","r")
labels = ["w_1", "w_2", "w_3", "w_4"]
for i in range(1,5):
    data = h5f[labels[i-1]].value
    dec = data[:,1]
    ra = data[:,0]
    ra_min, ra_max = ra.min(), ra.max()
    dec_min, dec_max = dec.min(), dec.max()
    pts_1 = [ra_min, ra_max, ra_max, ra_min, ra_min]
    pts_2 = [dec_min, dec_min, dec_max, dec_max, dec_min]

    img.axs[0][0].plot(pts_1, pts_2, label=labels[i-1])

h5f.close()
img.axs[0][0].legend()
img.save_img("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/pic/overlap_ori.png")
img.close()


# The catalog in grid
img = Image_Plot(fig_x=12, fig_y=9)
img.set_style()
img.plot_img(1, 1)

img.axs[0][0].scatter(ra_n, dec_n, s=3, alpha=0.5, label="CMASS_NORTH")
img.axs[0][0].scatter(ra_s, dec_s, s=3, alpha=0.5, label="CMASS_SOUTH")

h5f = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_grid.hdf5", "r")
labels = ["w_1", "w_2", "w_3", "w_4"]
for i in range(1, 5):
    ra = h5f["/background/w_%d/RA" % i].value
    dec = h5f["/background/w_%d/DEC" % i].value
    ra_bin = h5f["/background/w_%d/RA_bin" % i].value
    dec_bin = h5f["/background/w_%d/DEC_bin" % i].value
    shape = h5f["/background/w_%d" % i].attrs["grid_shape"]
    block_scale = h5f["/background/w_%d" % i].attrs["block_scale"]

    ra_min, ra_max = ra.min(), ra.max()
    dec_min, dec_max = dec.min(), dec.max()
    pts_1 = [ra_min, ra_max, ra_max, ra_min, ra_min]
    pts_2 = [dec_min, dec_min, dec_max, dec_max, dec_min]

    img.axs[0][0].plot(pts_1, pts_2, label=labels[i - 1])

h5f.close()
img.axs[0][0].legend()
img.save_img("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/pic/overlap.png")
img.close()



# The overlap
ras = [ra_s, ra_n]
decs = [dec_s, dec_n]
zs = [z_s,  z_n]
boss_label = ["CMASS_SOUTH", "CMASS_NORTH"]

img = Image_Plot(fig_x=8, fig_y=6)
img.set_style()
img.plot_img(2, 3)
h5f = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_grid.hdf5", "r")
sky_labels = ["w_1", "w_2", "w_3", "w_4"]

tag = [1, 3, 4]
data_tag = [0, 1, 0]
select_data = [[], []]
for j in range(2):
    for i in range(3):
        ra = h5f["/background/w_%d/RA" % tag[i]].value
        dec = h5f["/background/w_%d/DEC" % tag[i]].value
        ra_bin = h5f["/background/w_%d/RA_bin" % tag[i]].value
        dec_bin = h5f["/background/w_%d/DEC_bin" % tag[i]].value
        shape = h5f["/background/w_%d" % tag[i]].attrs["grid_shape"]
        block_scale = h5f["/background/w_%d" % tag[i]].attrs["block_scale"]

        ra_min, ra_max = ra.min(), ra.max()
        dec_min, dec_max = dec.min(), dec.max()

        pts_1 = [ra_min, ra_max, ra_max, ra_min, ra_min]
        pts_2 = [dec_min, dec_min, dec_max, dec_max, dec_min]

        img.axs[j][i].plot(pts_1, pts_2, linestyle="--",c="dimgrey", label=sky_labels[tag[i] - 1])

        idx1 = ras[data_tag[i]] >= ra_min
        idx2 = ras[data_tag[i]] <= ra_max
        idx3 = decs[data_tag[i]] >= dec_min
        idx4 = decs[data_tag[i]] <= dec_max
        idx = idx1 & idx2 & idx3 & idx4

        if j == 0:
            norm = plt.Normalize(vmin=numpy.min(zs[data_tag[i]][idx]), vmax=numpy.max(zs[data_tag[i]][idx]))
            cmap = plt.get_cmap('YlOrRd')
            img.axs[j][i].scatter(ras[data_tag[i]][idx], decs[data_tag[i]][idx], s=5, color=cmap(zs[data_tag[i]][idx]),
                                  label=boss_label[data_tag[i]])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm._A = []
            num = idx.sum()

            select_data[j].append([ras[data_tag[i]][idx], decs[data_tag[i]][idx], zs[data_tag[i]][idx]])
        else:
            idx_z1 = zs[data_tag[i]] >= 0.43
            idx_z2 = zs[data_tag[i]] <= 0.7

            idx_t = idx & idx_z1 & idx_z2

            norm = plt.Normalize(vmin=numpy.min(zs[data_tag[i]][idx_t]), vmax=numpy.max(zs[data_tag[i]][idx_t]))
            cmap = plt.get_cmap('YlOrRd')
            img.axs[j][i].scatter(ras[data_tag[i]][idx_t], decs[data_tag[i]][idx_t], s=5,
                                  color=cmap(zs[data_tag[i]][idx_t]), label=boss_label[data_tag[i]])
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm._A = []
            num = idx_t.sum()
            select_data[j].append([ras[data_tag[i]][idx_t], decs[data_tag[i]][idx_t], zs[data_tag[i]][idx_t]])

        img.tick_label(j, i, 1, "R.A.")
        img.tick_label(j, i, 0, "DEC.")
        img.axs[j][i].set_title(" Total BOSS galaxy: %d" % num)
        img.axs[j][i].legend()
        plt.colorbar(sm, ax=img.axs[j][i])

        # print(shape, block_scale)
        # print(dec_bin.min(), dec_bin.max(), dec_bin.shape)
        # print(ra_bin.min(), ra_bin.max(), ra_bin.shape)
        # print(dec_min, dec_max)
        # print(ra_min, ra_max, ra.shape, idx.sum())
        # print("\n")

img.subimg_adjust(0.3, 0.2)
img.save_img("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/pic/areas.pdf")
img.close()

h5f.close()


# The foreground
# write to the final hdf5 file for calculation
#  the CMASS catalog

h5f = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5","r")
redshift = h5f["/redshift"].value
distance = h5f["/distance"].value
h5f.close()

area_nm = ["w_1", "w_3", "w_4"]
h5f_cata = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_grid.hdf5", "a")
for i in range(3):
    ra_select = select_data[1][i][0]
    dec_select = select_data[1][i][1]
    dec_arc = select_data[1][i][1] / 180 * numpy.pi
    cos_dec_select = numpy.cos(dec_arc)
    z_select = select_data[1][i][2]
    dist_select = numpy.zeros_like(z_select)
    sp = [ra_select.shape[0], 0]

    for iz in range(sp[0]):
        tag = tool_box.find_near(redshift, z_select[iz])
        dist_select[iz] = distance[tag]

    h5f_cata["/foreground/%s/Z"] = z_select
    h5f_cata["/foreground/%s/Z"].attrs["shape"] = numpy.array(sp, dtype=numpy.intc)

    h5f_cata["/foreground/%s/DISTANCE"] = dist_select
    h5f_cata["/foreground/%s/DISTANCE"].attrs["shape"] = numpy.array(sp, dtype=numpy.intc)

    h5f_cata["/foreground/%s/RA"] = ra_select
    h5f_cata["/foreground/%s/RA"].attrs["shape"] = numpy.array(sp, dtype=numpy.intc)

    h5f_cata["/foreground/%s/DEC"] = dec_select
    h5f_cata["/foreground/%s/DEC"].attrs["shape"] = numpy.array(sp, dtype=numpy.intc)

    h5f_cata["/foreground/%s/COS_DEC"] = cos_dec_select
    h5f_cata["/foreground/%s/COS_DEC"].attrs["shape"] = numpy.array(sp, dtype=numpy.intc)
h5f_cata.close()