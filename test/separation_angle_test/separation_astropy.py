import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import time
import tool_box
import h5py
import numpy
from mpi4py import MPI
from astropy.coordinates import SkyCoord
from astropy import units as astro_unit


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
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/cata_result_ext_grid.hdf5"


h5f_fore = h5py.File(fore_path,"r")
# h5f_fore = h5py.File("/home/hkli/work/test/cluster_w_3.hdf5", "r")
RA_f = h5f_fore["/RA"].value
DEC_f = h5f_fore["/DEC"].value
COS_DEC_f = h5f_fore["/COS_DEC"].value
DISTANCE_f = h5f_fore["/DISTANCE"].value[:,0]
Z_f = h5f_fore["/Z"].value
h5f_fore.close()

h5f_data = h5py.File(data_path, "r")
RA_d = h5f_data["/%s/RA"%area_id].value[:,0]
DEC_d = h5f_data["/%s/DEC"%area_id].value[:,0]
COS_DEC_d = h5f_data["/%s/COS_DEC"%area_id].value[:,0]
DISTANCE_d = h5f_data["/%s/DISTANCE"%area_id].value[:,0]
Z_d = h5f_data["/%s/Z"%area_id].value[:,0]
radius_bin = h5f_data["/radius_bin"].value[:,0]
h5f_data.close()


galaxy_pos = SkyCoord(ra=RA_d*astro_unit.deg, dec=DEC_d*astro_unit.deg, frame='fk5')

fore_num = RA_f.shape[0]
gal_num = Z_d.shape[0]

background_label = numpy.arange(gal_num)

tasks = [i for i in range(fore_num)]
my_task = tool_box.allot(tasks, cpus)[rank]

# my_task = [int(argv[1])]

itemsize = MPI.INT.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = fore_num*itemsize
else:
    nbytes = 0
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
pair_num = numpy.ndarray(buffer=buf1, dtype=numpy.intc, shape=(fore_num,)) # array filled with zero


t1 = time.time()


for ig in my_task:

    my_RA_f = RA_f[ig]
    my_DEC_f = DEC_f[ig]
    my_COS_DEC_f = COS_DEC_f[ig]
    my_DISTANCE_f = DISTANCE_f[ig]
    my_Z_f = Z_f[ig]

    my_pos = SkyCoord(ra=my_RA_f * astro_unit.deg, dec=my_DEC_f * astro_unit.deg, frame='fk5')
    sep_angle = my_pos.separation(galaxy_pos).radian

    idx = Z_d >= my_Z_f + 0.1

    back_num = idx.sum()

    if back_num > 0:
        dist_tran = sep_angle[idx]*DISTANCE_d[idx]*1000*C_0_hat # h^{-1} Mpc
        idx_1 = dist_tran >= radius_bin[12]
        idx_2 = dist_tran < radius_bin[13]
        idx_t = idx_1 & idx_2

        # print(dist_tran[idx_t], dist_tran[idx_t]/0.7)
        #
        pair_num[ig] = idx_t.sum()
        #
        # pair_data = numpy.zeros((int(pair_num[ig]), 10))
        #
        # pair_data[:, 0] = RA_d[idx][idx_t]
        # pair_data[:, 1] = DEC_d[idx][idx_t]
        # pair_data[:, 2] = COS_DEC_d[idx][idx_t]
        # pair_data[:, 3] = DISTANCE_d[idx][idx_t]
        # pair_data[:, 4] = Z_d[idx][idx_t]
        # pair_data[:, 5] = background_label[idx][idx_t]
        #
        # pair_data[:, 6] = my_RA_f
        # pair_data[:, 7] = my_DEC_f
        # pair_data[:, 8] = my_DISTANCE_f
        # pair_data[:, 9] = my_Z_f

comm.Barrier()
t2 = time.time()

if rank == 0:

    # numpy.savez("/home/hkli/work/test/pair_data.npz", pair_data)

    numpy.savez("/home/hkli/work/test/ggl_mask.npz", pair_num)#,pair_data)

    fore_label = numpy.arange(fore_num)

    idx = pair_num > 0.01
    total_pair = pair_num[idx].sum()

    print(fore_num, idx.sum(), total_pair, pair_num[idx])
    print(fore_label[idx])

    h5f = h5py.File("/home/hkli/work/test/cluster_w_3_12.hdf5", "w")
    h5f["/RA"] = RA_f[idx]
    h5f["/DEC"] = DEC_f[idx]
    h5f["/COS_DEC"] = COS_DEC_f[idx]
    h5f["/DISTANCE"] = DISTANCE_f[idx]
    h5f["/Z"] = Z_f[idx]
    h5f.close()
    print(fore_num, gal_num, pair_num.sum(), t2-t1)
    print(radius_bin)

