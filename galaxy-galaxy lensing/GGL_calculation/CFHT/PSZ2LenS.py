import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import tool_box
import numpy
from mpi4py import MPI
import h5py
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units as astro_unit
import time



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


foreground_source = "PSZ2LenS"
area_id = "w_3"
h = 0.7
Om0 = 0.31
H0 = 100*h
C_0_hat = 2.99792458

# directories
parent_path = "/mnt/perc/hklee/CFHT/"
fg_cata_path = parent_path + "gg_lensing/data/foreground/%s/%s.hdf5"%(foreground_source, area_id)
bk_cata_path = parent_path + "catalog/cfht_cata/cfht_cata.hdf5"
result_path = parent_path + "result/%s/%s/"%(foreground_source, area_id)

# set up radius bins
radius_num = 24
radius_bin = tool_box.set_bin_log(0.1, 25.12, radius_num)

# read foreground
h5f = h5py.File(fg_cata_path,"r")
RA_fg = h5f["/3/RA"].value
DEC_fg = h5f["/3/DEC"].value
Z_fg = h5f["/3/Z"].value
DIST_fg = h5f["/3/DISTANCE"].value
DIST_INTEG_fg = h5f["/3/DISTANCE_INTEG"].value
h5f.close()

if rank == 0:
    print("Radius bin: ", radius_bin)
    print("fore Z: ", Z_fg)

# read background
ra_lb_c = 0
dec_lb_c = 1
flag_lb_c = 2
flux_rad_lb_c = 3
e1_lb_c = 4
e2_lb_c = 5
weight_lb_c = 6
fitclass_lb_c = 7
snr_lb_c = 8
mask_lb_c = 9
z_lb_c = 10
m_lb_c = 11
c_lb_c = 12
lp_mi_lb_c = 13
starflag_lb_c = 14
mag_lb_c = 15
z_min_lb_c = 16
z_max_lb_c = 17
odds_lb_c = 18

h5f = h5py.File(bk_cata_path,"r")
data_bk = h5f["/%s"%area_id].value
h5f.close()

# cut off
epsilon = 0.0000001

# weight > 0
idx_weight = data_bk[:, weight_lb_c] > 0 - epsilon

# # mask <= 1
# idx_mask = data_bk[:, mask_lb_c] <= 1 + epsilon
#
# # FITCLASS = 0
# idx_fitclass = numpy.abs(data_bk[:, fitclass_lb_c]) < 0 + epsilon
#
# # MAG < 24.7
# idx_mag = data_bk[:, mag_lb_c] < 24.7

# Redshift
idxz_1 = data_bk[:, z_lb_c] <= 1.2
idxz_2 = data_bk[:, z_lb_c] >= numpy.min(Z_fg)
idxz_3 = data_bk[:, odds_lb_c] > 0.8 - epsilon

cut_idx = idx_weight & idxz_1 & idxz_2 & idxz_3

RA_bk = data_bk[:,ra_lb_c][cut_idx]
DEC_bk = data_bk[:,dec_lb_c][cut_idx]
Z_bk = data_bk[:,z_lb_c][cut_idx]
Zmin_bk = data_bk[:,z_min_lb_c][cut_idx]
Zmax_bk = data_bk[:,z_max_lb_c][cut_idx]
ODDS_bk = data_bk[:,odds_lb_c][cut_idx]
num_bk = Z_bk.shape[0]

itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = num_bk*itemsize
else:
    nbytes = 0
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
DIST_INTEG_bk = numpy.ndarray(buffer=buf1, dtype='d', shape=(num_bk,))

m, n = divmod(num_bk, cpus)
my_pts_s, my_pts_e = m*rank, (rank+1)*m
if rank == cpus-1:
    my_pts_e += n
    print("Background: ",num_bk, m, n)

# for separation calculation
cosmos = FlatLambdaCDM(H0, Om0)

t1 = time.time()
for i in range(my_pts_s, my_pts_e):
    DIST_INTEG_bk[i] = cosmos.comoving_distance(Z_bk[i]).value*h/1000/C_0_hat
t2 = time.time()
comm.Barrier()

if rank == 0:

    print(DIST_INTEG_bk[-n:])
    print(t2-t1)


