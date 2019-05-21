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

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

t1 = time.time()

foreground_name = "CFHT_cluster"

area_id = int(argv[1])
# ra ~ foreground_ra +/- delta_ra
# input arcminute  -> degree
delta_ra = float(argv[2])/60
delta_dec = delta_ra
# grid nx, ny (even number)
nx = int(argv[3])
max_num = max(0,int(argv[4]))

ny = nx
grid_num = nx*ny

delta_z = 0.1
# degree, separation angle for shear estimation
radius = 10./60
# degree, smooth scale in the weight
smooth_len = 30./3600

total_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/mass_map/"
result_path = total_path + foreground_name + "/w_%d/"%area_id
if rank == 0:
    if not os.path.exists(total_path+foreground_name):
        os.mkdir(total_path+foreground_name)
    if not os.path.exists(result_path):
        os.mkdir(result_path)
comm.Barrier()

data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/"

fq = Fourier_Quad(10, 123)

# gamma_t, sig, gamma_x, sig, position angle, galaxy number
block_num = 6
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

#  read foreground clusters
h5f = h5py.File(data_path + "foreground/%s/w_%d_sub.hdf5"%(foreground_name, area_id), "r")
fore_ra = h5f["/RA"].value
fore_dec = h5f["/DEC"].value
fore_z = h5f["/Z"].value
fore_dist = h5f["/DISTANCE"].value
fore_cos_dec = h5f["/COS_DEC"].value
fore_n200 = h5f["/N200"].value
h5f.close()

fore_num = fore_ra.shape[0]

# read background catalog
names = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V", "COS_DEC", "DISTANCE"]

grid_id = [(i, j) for i in range(ny) for j in range(nx)]
my_grid = tool_box.allot(grid_id, cpus)[rank]

t2 = time.time()

for igal in range(min(max_num,fore_num)):

    it1 = time.time()

    h5f = h5py.File(data_path + "foreground/%s/w_%d_sub.hdf5" % (foreground_name, area_id), "r")

    redshift = h5f["/%d/Z" %igal].value
    data = numpy.zeros((redshift.shape[0], len(names)))
    data[:, 0] = redshift
    for i in range(1, len(names)):
        data_arr = h5f["/%d/%s" % (igal, names[i])].value
        if len(data_arr.shape) > 1:
            data[:, i] = data_arr[:, 0]
        else:
            data[:, i] = data_arr
    # the x-coord is opposite to RA
    data[:, names.index("G2")] = - data[:, names.index("G2")]
    h5f.close()

    result_path_ig = result_path + "source_%d/"%igal
    pic_path_ig = result_path_ig + "pic/"
    if rank == 0:
        if not os.path.exists(result_path_ig):
            os.mkdir(result_path_ig)

        if not os.path.exists(pic_path_ig):
            os.mkdir(pic_path_ig)

        for i in range(block_num * ny):
            for j in range(nx):
                result[i, j] = -99
    comm.Barrier()

    center = SkyCoord(ra=fore_ra[igal]*astro_unit.degree,dec=fore_dec[igal]*astro_unit.degree, frame='icrs')

    # set up bins for RA and DEC
    ra_bin = numpy.linspace(fore_ra[igal]-delta_ra, fore_ra[igal]+delta_ra, nx+1)
    dec_bin = numpy.linspace(fore_dec[igal]-delta_dec, fore_dec[igal]+delta_dec, ny+1)

    idx_z = data[:, names.index("Z")] >= fore_z[igal] + delta_z

    # sub-area
    mg_t = data[:, names.index("G1")][idx_z]
    mg_x = data[:, names.index("G2")][idx_z]
    mu = data[:, names.index("U")][idx_z]
    mn = data[:, names.index("N")][idx_z]
    ra = data[:, names.index("RA")][idx_z]
    dec = data[:, names.index("DEC")][idx_z]

    galaxy_pos = SkyCoord(ra=ra*astro_unit.degree, dec=dec*astro_unit.degree, frame='icrs')

    for grid in my_grid:
        iy, ix = grid

        grid_ra = (ra_bin[ix] + ra_bin[ix + 1]) / 2
        grid_dec = (dec_bin[iy] + dec_bin[iy + 1]) / 2
        grid_center = SkyCoord(ra=grid_ra*astro_unit.degree, dec=grid_dec*astro_unit.degree, frame='icrs')

        sep_angle = grid_center.separation(galaxy_pos).degree

        idx_sub = sep_angle <= radius

        weight = numpy.exp(-sep_angle[idx_sub]**2/2/smooth_len**2)

        mg_t_ = mg_t[idx_sub]*weight
        mg_x_ = mg_x[idx_sub]*weight
        mn_ = mn[idx_sub]*weight
        mu_ = mu[idx_sub]*weight
        mnu1 = mn_ + mu_
        mnu2 = mn_ - mu_

        result[iy + 4 * ny, ix] = idx_sub.sum()

        if idx_sub.sum() > 0:
            try:
                pic_nm = pic_path_ig + "%d_%d_t.png" % (iy, ix)
                tan_g, tan_g_sig = fq.fmin_g_new(g=mg_t_, nu=mnu1, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1, right=0.1,
                                                 fit_num=20)
                result[iy, ix] = tan_g
                result[iy + ny, ix] = tan_g_sig
            except:
                pass
            try:
                pic_nm = pic_path_ig + "%d_%d_x.png" % (iy, ix)
                cross_g, cross_g_sig = fq.fmin_g_new(g=mg_x_, nu=mnu2, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1,
                                                     right=0.1, fit_num=20)
                result[iy + 2 * ny, ix] = cross_g
                result[iy + 3 * ny, ix] = cross_g_sig
            except:
                pass

    if rank == 0:
        numpy.savez(result_path_ig + "result.npz", result, ra_bin, dec_bin,[fore_ra[igal], fore_dec[igal], fore_z[igal]])

        log_path = result_path_ig + "result.dat"

        paras = [["para", "area", str(area_id)], ["para", "foreground Z", str(fore_z[igal])],
                 ["para", "foreground RA", str(fore_ra[igal])],
                 ["para", "foreground Dec", str(fore_dec[igal])], ["para", "delta RA", str(delta_ra*60)],
                 ["para", "delta Dec", str(delta_dec*60)],
                 ["para", "ny", str(ny)], ["para", "nx", str(nx)], ["para", "delta Z", str(delta_z)],
                 ["para", "total galaxy", str(idx_z.sum())]]
        cmd = ["add" for i in range(len(paras))]
        tool_box.config(log_path, cmd, paras, write=True)
    it2 = time.time()
    print(rank, "%d galaxy: %.2f sec" %(igal, it2 - it1))
    comm.Barrier()
comm.Barrier()

t3 = time.time()
print(rank, " %.2f,%.2f"%(t2-t1, t3-t2))



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
