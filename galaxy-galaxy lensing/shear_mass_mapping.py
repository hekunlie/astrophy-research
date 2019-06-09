import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import time
from Fourier_Quad import Fourier_Quad
import tool_box
import plot_tool
import h5py
import numpy
from mpi4py import MPI
from astropy.coordinates import SkyCoord
from astropy import units as astro_unit
import matplotlib.pyplot as plt

# the RA & DEC are calculated in arcminute
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

t1 = time.time()

foreground_name = "CFHT_cluster"
deg2arcmin = 60
arcmin2rad = 1./60/180*numpy.pi

area_id = int(argv[1])
# ra ~ foreground_ra +/- delta_ra
# input arcminute  -> degree
delta_ra = float(argv[2])
delta_dec = delta_ra
# grid nx, ny (even number)
nx = int(argv[3])
max_num = max(0,int(argv[4]))

ny = nx
grid_num = nx*ny

delta_z = 0.1
# arcmin, separation angle for shear estimation
radius = 5
radius_sq = radius**2
# arcmin, smooth scale in the weight
smooth_len = 2

sgima_coeff = 388.283351

total_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/mass_map/"
result_path = total_path + foreground_name + "/w_%d/"%area_id
if rank == 0:
    if not os.path.exists(total_path+foreground_name):
        os.mkdir(total_path+foreground_name)
    if not os.path.exists(result_path):
        os.mkdir(result_path)
    if not os.path.exists(result_path + "pic/"):
        os.mkdir(result_path + "pic/")
comm.Barrier()

data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/"

fq = Fourier_Quad(10, 123)

# gamma_t, sig, gamma_x, sig, position angle, galaxy number
block_num = 10
itemsize = MPI.DOUBLE.Get_size()
element_num = nx*ny*block_num
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
result = numpy.ndarray(buffer=buf1, dtype='d', shape=(block_num*ny, nx))

comm.Barrier()
if rank == 0:
    if not os.path.exists(result_path):
        os.mkdir(result_path)
comm.Barrier()

# h5f = h5py.File(data_path + "redshift.hdf5","r")
# redshift_refer = h5f["/redshift"].value
# dist_refer = h5f["/distance"].value
# h5f.close()

#  read foreground clusters
h5f = h5py.File(data_path + "foreground/%s/w_%d_sub.hdf5"%(foreground_name, area_id), "r")
fore_dist = h5f["/DISTANCE"].value
fore_cos_dec = h5f["/COS_DEC"].value
fore_n200 = h5f["/N200"].value
fore_ra = h5f["/RA"].value*deg2arcmin
fore_dec = h5f["/DEC"].value*deg2arcmin
fore_z = h5f["/Z"].value
h5f.close()


fore_num = fore_ra.shape[0]

# read background catalog
names = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V", "COS_DEC", "DISTANCE"]

grid_id = [(i, j) for i in range(ny) for j in range(nx)]
my_grid = tool_box.allot(grid_id, cpus)[rank]

t2 = time.time()

for igal in range(min(max_num,fore_num)):

    it1 = time.time()
    # created by CFHT_cluster_cata.ipynb +/- 0.75 degree around each foreground point
    h5f = h5py.File(data_path + "foreground/%s/w_%d_sub.hdf5" % (foreground_name, area_id), "r")

    redshift_all = h5f["/%d/Z" %igal].value
    data = numpy.zeros((redshift_all.shape[0], len(names)))
    data[:, 0] = redshift_all
    for i in range(1, len(names)):
        data_arr = h5f["/%d/%s" % (igal, names[i])].value
        if len(data_arr.shape) > 1:
            data[:, i] = data_arr[:, 0]
        else:
            data[:, i] = data_arr
    # the x-coord is opposite to RA
    data[:, names.index("G2")] = - data[:, names.index("G2")]
    h5f.close()

    # set up redshift bin
    redshift_bin = [fore_z[igal] + delta_z]
    while redshift_bin[-1] <= 0:
        redshift_bin.append(redshift_bin[-1]+0.3)
    redshift_bin.append(20)

    result_path_ig = result_path + "source_%d/"%igal
    pic_path_ = result_path + "pic/source_%d/"%igal
    pic_path_ig = result_path + "pic/source_%d/pic/"%igal
    if rank == 0:
        if not os.path.exists(result_path_ig):
            os.mkdir(result_path_ig)
        if not os.path.exists(pic_path_):
            os.mkdir(pic_path_)
        if not os.path.exists(pic_path_ig):
            os.mkdir(pic_path_ig)

        total_result = numpy.zeros(((len(redshift_bin) - 1) * block_num * ny, nx))

    # center = SkyCoord(ra=fore_ra[igal]*astro_unit.degree,dec=fore_dec[igal]*astro_unit.degree, frame='fk5')

    # set up bins for RA and DEC
    ra_bin = numpy.linspace(fore_ra[igal]-delta_ra, fore_ra[igal]+delta_ra, nx+1)
    dec_bin = numpy.linspace(fore_dec[igal]-delta_dec, fore_dec[igal]+delta_dec, ny+1)
    ra_bin_at_dec = ra_bin*numpy.cos(dec_bin*arcmin2rad)

    idx_z = data[:, names.index("Z")] >= fore_z[igal] + delta_z

    # sub-area
    redshift = data[:, names.index("Z")][idx_z]
    dist = data[:,names.index("DISTANCE")][idx_z]

    sigma_crit = sgima_coeff*dist/(dist - fore_dist[igal])/fore_dist[igal]*(1+fore_z[igal])

    ra = data[:, names.index("RA")][idx_z]*deg2arcmin
    dec = data[:, names.index("DEC")][idx_z]*deg2arcmin
    cos_dec = data[:,names.index("COS_DEC")][idx_z]
    ra_at_dec = ra*cos_dec

    mg_t = data[:, names.index("G1")][idx_z]
    mg_x = data[:, names.index("G2")][idx_z]
    mu = data[:, names.index("U")][idx_z]
    mn = data[:, names.index("N")][idx_z]

    galaxy_pos = SkyCoord(ra=ra * astro_unit.arcmin, dec=dec * astro_unit.arcmin, frame='fk5')

    ra_min, ra_max = ra_bin.min(), ra_bin.max()
    dec_min, dec_max = dec_bin.min(), dec_bin.max()

    idx_ra_s1 = ra >= ra_min - radius
    idx_ra_s2 = ra <= ra_max + radius
    idx_dec_s1 = dec >= dec_min - radius
    idx_dec_s2 = dec <= dec_max + radius

    idx_p = idx_ra_s1 & idx_ra_s2 & idx_dec_s1 & idx_dec_s2

    for ir in range(len(redshift_bin) - 1):

        idx_z_b1 = redshift >= redshift_bin[ir]
        idx_z_b2 = redshift < redshift_bin[ir+1]

        idx_z_sub = idx_z_b1 & idx_z_b2

        idx_s = idx_p & idx_z_sub

        if rank == 0:

            img = plot_tool.Image_Plot(fig_x=12, fig_y=12)
            img.create_subfig(1, 2)

            img.axs[0][0].scatter(fore_ra[igal], fore_dec[igal], s=200, facecolors="none", edgecolors="r", marker="*")
            for i in range(ny + 1):
                img.axs[0][0].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--",
                                   alpha=0.5, linewidth=0.3)
            for j in range(nx + 1):
                img.axs[0][0].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black", linestyle="--",
                                   alpha=0.5, linewidth=0.3)
            if idx_s.sum() > 0:
                img.axs[0][0].scatter(ra[idx_s], dec[idx_s],s=3)
                img.axs[0][1].hist(redshift[idx_s], 100)

            img.tick_label(0,0,0,"DEC (arcmin)")
            img.tick_label(0,0,1,"RA (arcmin)")
            img.axs[0][0].set_title("%d galaxies (%.2f ~ %.2f)"%(idx_s.sum(),redshift_bin[ir], redshift_bin[ir+1]))

            img.tick_label(0,1,0,"Number")
            img.tick_label(0,1,1,"Z")
            img.save_img(result_path_ig+"s%d_density_%d.png"%(igal,ir))
            img.close_img()

            for i in range(block_num * ny):
                for j in range(nx):
                    result[i, j] = -2

        comm.Barrier()

        for grid in my_grid:
            iy, ix = grid

            grid_ra = (ra_bin[ix] + ra_bin[ix + 1]) / 2
            grid_dec = (dec_bin[iy] + dec_bin[iy + 1]) / 2
            grid_center = SkyCoord(ra=grid_ra*astro_unit.arcmin, dec=grid_dec*astro_unit.arcmin, frame='fk5')
            sep_angle = grid_center.separation(galaxy_pos).arcmin

            idx_sep = sep_angle <= radius
            idx_sub = idx_sep & idx_z_sub
            source_num = idx_sub.sum()

            weight = numpy.exp(-sep_angle[idx_sub]**2/2/smooth_len ** 2) * 40
            #
            # grid_ra = (ra_bin_at_dec[ix] + ra_bin_at_dec[ix + 1]) / 2
            # grid_dec = (dec_bin[iy] + dec_bin[iy + 1]) / 2
            # sep_angle = (ra_at_dec - grid_ra)**2 + (dec - grid_dec)**2
            # idx_sub = sep_angle <= radius_sq
            # weight = numpy.exp(-sep_angle[idx_sub]/2/smooth_len**2)*40
            sigma_crit_ = sigma_crit[idx_sub]
            mg_t_ = mg_t[idx_sub]*weight
            mg_x_ = mg_x[idx_sub]*weight
            mn_ = mn[idx_sub]*weight
            mu_ = mu[idx_sub]*weight
            mnu1 = mn_ + mu_
            mnu2 = mn_ - mu_
            # numpy.savez(total_path+"temp/%d_%d_%d_%d.npz"%(igal,iy,ix,ir), mg_t_, mg_x_, mn_, mu_, mnu1, mnu2)
            result[iy + (block_num - 1) * ny, ix] = source_num

            if source_num > 0:
                # try:
                #     pic_nm = pic_path_ig + "%d_%d_t_%d.png" % (iy, ix, ir)
                #     tan_g, tan_g_sig,coeff1 = fq.fmin_g_new(g=mg_t_, nu=mnu1, bin_num=8, scale=100, pic_path=pic_nm, left=-150,
                #                                      right=150, fit_num=20,chi_gap=100)
                # except:
                #     print("%d, %d, [%d, %d], Bad fitting g1" % (rank, igal, iy, ix))
                # try:
                #     pic_nm = pic_path_ig + "%d_%d_x_%d.png" % (iy, ix, ir)
                #     cross_g, cross_g_sig,coeff2 = fq.fmin_g_new(g=mg_x_, nu=mnu2, bin_num=8, scale=100, pic_path=pic_nm, left=-150,
                #                                          right=150, fit_num=20,chi_gap=100)
                # except:
                #     print("%d, %d, [%d, %d], Bad fitting g2" % (rank, igal, iy, ix))
                fig = plt.figure(figsize=(6,4.8))
                ax1 = fig.add_subplot(221)
                ax2 = fig.add_subplot(222)
                ax3 = fig.add_subplot(223)
                ax4 = fig.add_subplot(224)

                tan_g, tan_g_sig = fq.fmin_g_new(g=mg_t_*sigma_crit_, nu=mnu1, bin_num=8, scale=100,
                                                 fig_ax=ax1, left=-1000,right=1000, fit_num=20, chi_gap=30)[:2]
                result[iy, ix] = tan_g
                result[iy + ny, ix] = tan_g_sig

                cross_g, cross_g_sig = fq.fmin_g_new(g=mg_x_*sigma_crit_, nu=mnu2, bin_num=8, scale=100,
                                                     fig_ax=ax2, left=-1000,right=1000, fit_num=20, chi_gap=30)[:2]
                result[iy + 2*ny, ix] = cross_g
                result[iy + 3*ny, ix] = cross_g_sig

                result[iy + 4*ny, ix] = sigma_crit_.mean()

                # the tangential shear
                tan_g, tan_g_sig = fq.fmin_g_new(g=mg_t_, nu=mnu1, bin_num=8, scale=100,
                                                 fig_ax=ax3, left=-0.2,right=0.2, fit_num=20, chi_gap=30)[:2]
                result[iy + 5*ny, ix] = tan_g
                result[iy + 6*ny, ix] = tan_g_sig

                cross_g, cross_g_sig = fq.fmin_g_new(g=mg_x_, nu=mnu2, bin_num=8, scale=100,
                                                     fig_ax=ax4, left=-0.2,right=0.2, fit_num=20, chi_gap=30)[:2]
                result[iy + 7*ny, ix] = cross_g
                result[iy + 8*ny, ix] = cross_g_sig

                pic_nm = pic_path_ig + "%d_%d_%d.png" % (iy, ix, ir)
                plt.savefig(pic_nm, bbox_inches='tight')
                plt.close()

            else:
                print("%d, EMPTY %d, [%d, %d], %f"%(rank, igal, iy, ix,radius))
        comm.Barrier()
        if rank == 0:
            total_result[ir*block_num * ny:(ir+1)*block_num * ny,0:nx] = result
    comm.Barrier()
    if rank == 0:
        numpy.savez(result_path_ig + "result.npz", total_result, ra_bin, dec_bin, redshift_bin,[fore_ra[igal], fore_dec[igal], fore_z[igal]])

        log_path = result_path_ig + "result.dat"
        if os.path.exists(log_path):
            os.remove(log_path)

        paras = [["para", "area", str(area_id)], ["para", "foreground Z", str(fore_z[igal])],
                 ["para", "foreground RA (arcmin)", str(fore_ra[igal])], ["para", "foreground Dec (arcmin)", str(fore_dec[igal])],
                 ["para", "RA width (arcmin)", str(2*delta_ra)],
                 ["para", "Dec width (arcmin)", str(2*delta_dec)],
                 ["para", "delta RA (arcmin)", str(ra_bin[1]-ra_bin[0])],
                 ["para", "delta Dec (arcmin)", str(dec_bin[1] - dec_bin[0])],
                 ["para", "max radius (arcsec)", str(radius*60)], ["para", "smooth radius (arcsec)", str(smooth_len*60)],
                 ["para", "ny", str(ny)], ["para", "nx", str(nx)], ["para", "delta Z", str(delta_z)],["para","block num",str(block_num)]]
        cmd = ["add" for i in range(len(paras))]
        tool_box.config(log_path, cmd, paras, write=True)
    it2 = time.time()
    # print(rank, "%d galaxy: %.2f sec" %(igal, it2 - it1))
    comm.Barrier()
comm.Barrier()

t3 = time.time()
print(rank, "total: %.2f sec,%.2f sec"%(t2-t1, t3-t2))



# # rotation
# center = SkyCoord(ra=foreground_ra*astro_unit.degree,
#                   dec=foreground_dec*astro_unit.degree, frame='icrs')
#
# background_pos = SkyCoord(ra=data[:, names.index("RA")]*astro_unit.degree,
#                           dec=data[:, names.index("DEC")]*astro_unit.degree, frame='icrs')
#
# position_angle = center.position_angle(background_pos).radian
#
# # diff_ra = data[:, names.index("RA")] - foreground_ra
# # diff_dec = data[:, names.index("DEC")] - foreground_dec
#
# # sin_theta = numpy.sin(position_angle)
# # cos_theta = numpy.cos(position_angle)
# #
# # sin_2theta = numpy.sin(2*position_angle)
# # cos_2theta = numpy.cos(2*position_angle)
# #
# # sin_4theta = numpy.sin(4*position_angle)
# # cos_4theta = numpy.cos(4*position_angle)
#
# #sigma_crit = data[:, names.index("DISTANCE")]/(data[:, names.index("DISTANCE")] - foreground_dist)/foreground_dist*(1+foreground_z)*388.283351
# sigma_crit = 1
# # mg_t = (data[:, names.index("G1")]*cos_2theta - data[:, names.index("G2")]*sin_2theta)*sigma_crit
# # mg_x = (data[:, names.index("G1")]*sin_2theta + data[:, names.index("G2")]*cos_2theta)*sigma_crit
# # mu = (data[:, names.index("U")]*cos_4theta - data[:, names.index("V")]*sin_4theta)*sigma_crit
#
# mg_t = data[:, names.index("G1")]*sigma_crit
# mg_x = data[:, names.index("G2")]*sigma_crit
# mu = data[:, names.index("U")]*sigma_crit
#
# mn = data[:, names.index("N")]*sigma_crit
#
# t3 = time.time()
#
# for grid in my_grid:
#     iy, ix = grid
#
#     grid_pos = SkyCoord(ra=(ra_bin[ix] + ra_bin[ix+1])/2 * astro_unit.degree,
#                         dec=(dec_bin[iy] + dec_bin[iy+1])/2 * astro_unit.degree, frame='icrs')
#
#     # radian
#     sep_angle = center.position_angle(grid_pos).radian
#
#     idx_1 = data[:, names.index("RA")] >= ra_bin[ix]
#     idx_2 = data[:, names.index("RA")] < ra_bin[ix+1]
#     idx_3 = data[:, names.index("DEC")] >= dec_bin[iy]
#     idx_4 = data[:, names.index("DEC")] < dec_bin[iy+1]
#
#     idx_5 = data[:, names.index("Z")] >= foreground_z + delta_z
#
#     idx = idx_1 & idx_2 & idx_3 & idx_4
#     idx_z = idx_5 & idx
#
#     sub_num = idx.sum()
#     sub_num_back = idx_z.sum()
#
#     mg_t_ = mg_t[idx]
#     mg_x_ = mg_x[idx]
#     mn_ = mn[idx]
#     mu_ = mu[idx]
#     mnu1 = mn_ + mu_
#     mnu2 = mn_ - mu_
#
#     pic_nm = result_path + "mass_map/pic/%d_%d_t.png"%(iy,ix)
#     tan_g, tan_g_sig = fq.fmin_g_new(g=mg_t_, nu=mnu1, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1, right=0.1, fit_num=20)
#
#     pic_nm = result_path + "mass_map/pic/%d_%d_x.png"%(iy,ix)
#     cross_g, cross_g_sig = fq.fmin_g_new(g=mg_x_, nu=mnu2, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1, right=0.1, fit_num=20)
#
#     mg_t_ = mg_t[idx_z]
#     mg_x_ = mg_x[idx_z]
#     mn_ = mn[idx_z]
#     mu_ = mu[idx_z]
#     mnu1 = mn_ + mu_
#     mnu2 = mn_ - mu_
#
#     pic_nm = result_path + "mass_map/pic/%d_%d_t_back.png"%(iy,ix)
#     tan_g_back, tan_g_sig_back = fq.fmin_g_new(g=mg_t_, nu=mnu1, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1, right=0.1, fit_num=20)
#
#     pic_nm = result_path + "mass_map/pic/%d_%d_x_back.png"%(iy,ix)
#     cross_g_back, cross_g_sig_back = fq.fmin_g_new(g=mg_x_, nu=mnu2, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1, right=0.1, fit_num=20)
#
#     result[iy, ix] = tan_g_back
#     result[iy + ny, ix] = tan_g_sig_back
#     result[iy + 2*ny, ix] = cross_g_back
#     result[iy + 3*ny, ix] = cross_g_sig_back
#
#     result[iy + 4*ny, ix] = sub_num_back
#
#     result[iy + 5*ny, ix] = tan_g
#     result[iy + 6*ny, ix] = tan_g_sig
#     result[iy + 7*ny, ix] = cross_g
#     result[iy + 8*ny, ix] = cross_g_sig
#
#     result[iy + 9*ny, ix] = sub_num
#
#     result[iy + 10*ny, ix] = sep_angle
#
# comm.Barrier()
# t4 = time.time()
# print(rank, " %.2f,%.2f,%.2f"%(t2-t1, t3-t2, t4-t3))
#
# if rank == 0:
#     numpy.savez(result_path + "mass_map/result.npz", result, ra_bin, dec_bin,[foreground_ra, foreground_dec, foreground_z])
#
#     log_path = result_path + "mass_map/result.dat"
#
#     paras = [["para", "area", str(area_id)], ["para", "foreground Z", str(foreground_z)], ["para", "foreground RA", str(foreground_ra)],
#              ["para", "foreground Dec", str(foreground_dec)], ["para", "delta RA", str(delta_ra)], ["para", "delta Dec", str(delta_dec)],
#              ["para", "ny", str(ny)], ["para", "nx", str(nx)], ["para", "delta Z", str(delta_z)]]
#     cmd = ["add" for i in range(len(paras))]
#     tool_box.config(log_path, cmd, paras, write=True)
