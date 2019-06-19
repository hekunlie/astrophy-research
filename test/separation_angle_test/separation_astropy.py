import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import time
import tool_box
import plot_tool
import h5py
import numpy
from mpi4py import MPI
from astropy.coordinates import SkyCoord
from astropy import units as astro_unit
import matplotlib.pyplot as plt


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

deg2arcmin = 60
deg2rad = 1./180*numpy.pi

C_0_hat = 2.99792458
coeff = 0.18 / C_0_hat / numpy.pi
coeff_inv = C_0_hat * numpy.pi / 0.18

area_id = "w_3"
fore_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/foreground/CFHT_cluster/%s.hdf5"%area_id
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_cut.hdf5"

h5f_fore = h5py.File(fore_path,"r")
RA_f = h5f_fore["/RA"].value
DEC_f = h5f_fore["/DEC"].value
COS_DEC_f = h5f_fore["/COS_DEC"].value
DISTANCE_f = h5f_fore["/DISTANCE"].value[:,0]
Z_f = h5f_fore["/Z"].value
h5f_fore.close()

h5f_data = h5py.File(data_path, "r")
RA_d = h5f_data["/%s/RA"%area_id].value
DEC_d = h5f_data["/%s/DEC"%area_id].value
COS_DEC_d = h5f_data["/%s/COS_DEC"%area_id].value
DISTANCE_d = h5f_data["/%s/DISTANCE"%area_id].value[:,0]
Z_d = h5f_data["/%s/Z"%area_id].value
h5f_data.close()

galaxy_pos = SkyCoord(ra=RA_d*astro_unit.deg, dec=DEC_d*astro_unit.deg, frame='fk5')

fore_num = RA_f.shape[0]
tasks = [i for i in range(fore_num)]
my_task = tool_box.allot(tasks, cpus)[rank]

for ig in my_task:

    my_RA_f = RA_f[ig]
    my_DEC_f = DEC_f[ig]
    my_COS_DEC_f = COS_DEC_f[ig]
    my_DISTANCE_f = DISTANCE_f[ig]
    my_Z_f = Z_f[ig]

    my_pos = SkyCoord(ra=my_RA_f * astro_unit.deg, dec=my_DEC_f * astro_unit.deg, frame='fk5')
    sep_angle = my_pos.separation(galaxy_pos).radian

    idx = Z_d > my_Z_f + 0.1

    back_num = idx.sum()

    if back_num > 0:
        dist_tran = sep_angle[idx]*DISTANCE_d[idx]*1000*3 # h^{-1} Mpc
        idx_1 = dist_tran >= 0.04
        idx_2 = dist_tran < 0.0631048
        idx_t = idx_1 & idx_2


