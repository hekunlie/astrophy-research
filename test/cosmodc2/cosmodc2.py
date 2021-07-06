from sys import path
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py
from mpi4py import MPI
# import FQlib
import time
# import c4py
from astropy.cosmology import FlatLambdaCDM
from astropy import units
from astropy.coordinates import SkyCoord


# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()
# numprocs = comm.Get_size()
numprocs = 1# comm.Get_size()
rank = 0


h0 = 0.71
H0 = 100*h0
omg_cm0 = 0.2648
omg_bm0 = 0.0448
# omg_m0 = omg_cm0 + omg_bm0
omg_m0 = 0.2648

cosmos = FlatLambdaCDM(H0,omg_m0)

radius_bin_num = 13
radius_bin = tool_box.set_bin_log(0.1, 30, radius_bin_num+1)
if rank == 0:
    print("Radius bin: ", radius_bin)

h5f = h5py.File("/home/hklee/work/cosmoDC2/c_data.h5","r")
halo_dec = h5f["/dec"][()]
halo_ra = h5f["/ra"][()]
halo_mass = h5f["/halo_mass"][()]
halo_z = h5f["/redshift"][()]
# halo_id = h5f["/halo_id"][()].astype(dtype=numpy.intc)
h5f.close()

h5f = h5py.File("/home/hklee/work/cosmoDC2/g_data.h5","r")
ori_src_dec = h5f["/dec"][()]
ori_src_ra = h5f["/ra"][()]
ori_src_z = h5f["/redshift"][()]
ori_src_g1 = h5f["/shear_1"][()]
ori_src_g2 = h5f["/shear_2"][()]
# ori_src_halo_id = h5f["/halo_id"][()].astype(dtype=numpy.intc)
h5f.close()

idx1 = numpy.abs(ori_src_g1) <= 0.05
idx2 = numpy.abs(ori_src_g2) <= 0.05
idx = idx1 & idx2

src_dec = ori_src_dec[idx]
src_ra = ori_src_ra[idx]
src_z = ori_src_z[idx]
src_g1 = ori_src_g1[idx]
src_g2 = ori_src_g2[idx]
# src_halo_id = ori_src_halo_id[idx]
print(idx.sum(), idx.shape[0])

com_dist_halo = cosmos.comoving_distance(halo_z).value*h0
com_dist_src = cosmos.comoving_distance(src_z).value*h0

src_position = SkyCoord(ra = src_ra*units.deg, dec = src_dec*units.deg,frame="fk5")
# halo_position = SkyCoord(ra=halo_ra*units.deg, dec=halo_dec*units.deg)

halo_num = halo_dec.shape[0]
my_halo = tool_box.alloc([i for i in range(halo_num)], numprocs)[rank]


my_cache = numpy.zeros((4, radius_bin_num))
#
# itemsize = MPI.DOUBLE.Get_size()
#
# if rank == 0:
#     # bytes for 10 double elements
#     nbytes = 4*radius_bin_num*itemsize
# else:
#     nbytes = 0
#
# win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
# buf1, itemsize = win1.Shared_query(0)
# total_cache = numpy.ndarray(buffer=buf1, dtype='d', shape=(4, radius_bin_num))


for hl_id in my_halo:
    t1 = time.time()
    halo_pos_i = SkyCoord(ra=halo_ra[hl_id]*units.deg, dec=halo_dec[hl_id]*units.deg,frame="fk5")
    idxz = src_z >= halo_z[hl_id] + 0.1

    sep_ang = halo_pos_i.separation(src_position).radian
    sep_radius = sep_ang*com_dist_halo[hl_id]

    idx1 = sep_radius <= radius_bin[-1]

    idx = idx1 & idxz

    sep_radius_need = sep_radius[idx]
    com_dist_src_need = com_dist_src[idx]
    src_z_need = src_z[idx]

    src_g1_need = src_g1[idx]
    src_g2_need = src_g2[idx]

    pos_ang_need = halo_pos_i.position_angle(src_position[idx]).radian
    cos_2theta = numpy.cos(2 * pos_ang_need)
    sin_2theta = numpy.sin(2 * pos_ang_need)


    cri_density = com_dist_src_need/com_dist_halo[hl_id]/(com_dist_src_need - com_dist_halo[hl_id])/(1+src_z_need)*1662895.2007121066

    gt = (cos_2theta*src_g1_need - sin_2theta*src_g2_need)#*cri_density
    gx = (sin_2theta*src_g1_need + cos_2theta*src_g2_need)#*cri_density
    print(cri_density)
    for ir in range(radius_bin_num):

        idx1 = sep_radius_need >= radius_bin[ir]
        idx2 = sep_radius_need < radius_bin[ir+1]
        idx = idx1 & idx2

        my_cache[0,ir] += numpy.sum(gt[idx])
        my_cache[1,ir] += numpy.sum(gx[idx])
        my_cache[2,ir] += numpy.sum(sep_radius_need[idx])
        my_cache[3,ir] += idx.sum()
        print(hl_id, ir, idx.sum())
    t2 = time.time()
    print(hl_id, t2-t1)
    delta_sigma_t = my_cache[0]/my_cache[3]
    delta_sigma_x = my_cache[1]/my_cache[3]
    radius = my_cache[2]/my_cache[3]
    print(delta_sigma_t)
    print(delta_sigma_x)
    print(radius)
    print(radius_bin)

# comm.Barrier()
# #
# # for i in range(numprocs):
# #     if i == rank:
# #         total_cache += my_cache
# #     comm.Barrier()
# # comm.Barrier()
#
# if rank == 0:
#     # numpy.savez("./cache.npz",total_cache)
#     total_cache = numpy.load("./cache.npz")["arr_0"]
#
#     delta_sigma_t = total_cache[0]/total_cache[3]
#     delta_sigma_x = total_cache[1]/total_cache[3]
#     radius = total_cache[2]/total_cache[3]
#
#     img = Image_Plot()
#     img.subplots(1,2)
#     img.axs[0][0].plot(radius, delta_sigma_t)
#     img.axs[0][1].plot(radius, delta_sigma_x)
#     for i in range(2):
#         img.set_label(0,i,1, "Radius Mpc/h")
#
#         img.axs[0][i].set_xscale("log")
#         img.axs[0][i].set_yscale("log")
#     img.save_img("./result.pdf")
# comm.Barrier()





