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

result_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/"

data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/"

area_id = int(argv[1])
foreground_z = float(argv[2])
# degree
foreground_ra = float(argv[3])
foreground_dec = float(argv[4])
# ra ~ foreground_ra +/- delta_ra
# input arcminute  -> degree
delta_ra = float(argv[5])/60
delta_dec = delta_ra
# grid nx, ny (even number)
nx = int(argv[6])
ny = nx
grid_num = nx*ny

delta_z = 0.1

h5f = h5py.File(data_path + "redshift.hdf5", "r")
redshift_refer = h5f["/redshift"].value
dist_refer = h5f["/distance"].value
h5f.close()

tag = tool_box.find_near(redshift_refer, foreground_z)
foreground_dist = dist_refer[tag]

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
    for i in range(block_num*ny):
        for j in range(nx):
            result[i,j] = 0
comm.Barrier()

ra_bin = numpy.linspace(foreground_ra-delta_ra, foreground_ra+delta_ra, nx+1)
dec_bin = numpy.linspace(foreground_dec-delta_dec, foreground_dec+delta_dec, ny+1)


# selection
h5f = h5py.File(data_path + "cata_result_ext_cut.hdf5", "r")

names = ["Z", "RA", "DEC", "G1", "G2", "N", "U", "V", "COS_DEC", "DISTANCE"]

z = h5f["/w_%d/Z" % area_id].value
idx_z = z >= foreground_z + delta_z

ra = h5f["/w_%d/RA" % area_id].value
idx_ra_1 = ra >= ra_bin[0]
idx_ra_2 = ra < ra_bin[-1]

dec = h5f["/w_%d/DEC" % area_id].value
idx_dec_1 = dec >= dec_bin[0]
idx_dec_2 = dec < dec_bin[-1]

idx = idx_z & idx_ra_1 & idx_ra_2 & idx_dec_1 & idx_dec_2

background_num = idx.sum()
data = numpy.zeros((background_num, len(names)))
for tag, nm in enumerate(names):
    data_arr = h5f["/w_%d/%s" % (area_id, nm)].value[idx]
    if len(data_arr.shape) > 1:
        data[:, tag] = data_arr[:, 0]
    else:
        data[:, tag] = data_arr
# the x-coord is opposite to RA
data[:, names.index("G2")] = - data[:, names.index("G2")]
h5f.close()


grid_id = [(i, j) for i in range(ny) for j in range(nx)]
my_grid = tool_box.allot(grid_id, cpus)[rank]

t2 = time.time()

# rotation
center = SkyCoord(ra=foreground_ra*astro_unit.degree,
                  dec=foreground_dec*astro_unit.degree, frame='icrs')

background_pos = SkyCoord(ra=data[:, names.index("RA")]*astro_unit.degree,
                          dec=data[:, names.index("DEC")]*astro_unit.degree, frame='icrs')

sep_angle = center.position_angle(background_pos).radian

# diff_ra = data[:, names.index("RA")] - foreground_ra
# diff_dec = data[:, names.index("DEC")] - foreground_dec

sin_theta = numpy.sin(sep_angle)
cos_theta = numpy.cos(sep_angle)

sin_2theta = 2*sin_theta*cos_theta
cos_2theta = cos_theta**2 - sin_theta**2

sin_4theta = 2*sin_2theta*cos_2theta
cos_4theta = cos_2theta**2 - sin_2theta**2

sigma_crit = data[:, names.index("DISTANCE")]/(data[:, names.index("DISTANCE")] - foreground_dist)/foreground_dist*(1+foreground_z)*388.283351

mg_t = (data[:, names.index("G1")]*cos_2theta - data[:, names.index("G2")]*sin_2theta)*sigma_crit
mg_x = (data[:, names.index("G1")]*sin_2theta + data[:, names.index("G2")]*cos_2theta)*sigma_crit
mu = (data[:, names.index("U")]*cos_4theta - data[:, names.index("V")]*sin_4theta)*sigma_crit
mn = data[:, names.index("N")]*sigma_crit

t3 = time.time()

for grid in my_grid:
    iy, ix = grid

    grid_pos = SkyCoord(ra=(ra_bin[ix] + ra_bin[ix+1])/2 * astro_unit.degree,
                              dec=(dec_bin[iy] + dec_bin[iy+1])/2 * astro_unit.degree, frame='icrs')

    # radian
    sep_angle = center.position_angle(grid_pos).radian

    idx_1 = data[:, names.index("RA")] >= ra_bin[ix]
    idx_2 = data[:, names.index("RA")] < ra_bin[ix+1]
    idx_3 = data[:, names.index("DEC")] >= dec_bin[iy]
    idx_4 = data[:, names.index("DEC")] < dec_bin[iy+1]
    idx = idx_1 & idx_2 & idx_3 & idx_4
    sub_num = idx.sum()

    mg_t_ = mg_t[idx]
    mg_x_ = mg_x[idx]
    mn_ = mn[idx]
    mu_ = mu[idx]
    mnu1 = mn_ + mu_
    mnu2 = mn_ - mu_

    pic_nm = result_path + "mass_map/pic/%d_%d_t.png"%(iy,ix)
    tan_g, tan_g_sig = fq.fmin_g_new(g=mg_t_, nu=mnu1, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1, right=0.1, fit_num=20)

    pic_nm = result_path + "mass_map/pic/%d_%d_x.png"%(iy,ix)
    cross_g, cross_g_sig = fq.fmin_g_new(g=mg_x_, nu=mnu2, bin_num=10, scale=100, pic_path=pic_nm, left=-0.1, right=0.1, fit_num=20)

    result[iy, ix] = tan_g
    result[iy + ny, ix] = tan_g_sig
    result[iy + 2*ny, ix] = cross_g
    result[iy + 3*ny, ix] = cross_g_sig
    result[iy + 4*ny, ix] = sub_num
    result[iy + 5*ny, ix] = sep_angle

comm.Barrier()
t4 = time.time()
print(rank, " %.2f,%.2f,%.2f"%(t2-t1, t3-t2, t4-t3))

if rank == 0:
    numpy.savez(result_path + "mass_map/result.npz", result, ra_bin, dec_bin,[foreground_ra, foreground_dec, foreground_z])

    log_path = result_path + "mass_map/result.dat"

    paras = [["para", "area", str(area_id)], ["para", "foreground Z", str(foreground_z)], ["para", "foreground RA", str(foreground_ra)],
             ["para", "foreground Dec", str(foreground_dec)], ["para", "delta RA", str(delta_ra)], ["para", "delta Dec", str(delta_dec)],
             ["para", "ny", str(ny)], ["para", "nx", str(nx)], ["para", "delta Z", str(delta_z)]]
    cmd = ["add" for i in range(len(paras))]
    tool_box.config(log_path, cmd, paras, write=True)
