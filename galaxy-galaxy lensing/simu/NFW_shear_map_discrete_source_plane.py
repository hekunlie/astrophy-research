from sys import path
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import GGLensing_tool
import numpy
import h5py
import galsim
from Fourier_Quad import Fourier_Quad
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


fq = Fourier_Quad(12, 123)

# cosmology
omega_m0 = 0.31
omega_lam0 = 1 - omega_m0
h = 0.6735
C_0_hat = 2.99792458
H_0 = 100 * h
coeff = 1000 * C_0_hat / h

coeff_crit = C_0_hat ** 2 / 4 / numpy.pi / 6.674

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)

# Halo parameters
Mass = 10 ** 14  # M_sun/h
conc = 3.5  # concentration
len_z = 0.3  # redshift
halo_position = galsim.PositionD(0, 0)  # arcsec
com_dist_len = cosmos.comoving_distance(len_z).value * h  # Mpc/h
print("Lens plane at z = %.2f, %.5f Mpc/h" % (len_z, com_dist_len))

# lens profile
nfw = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)

# source parameters
src_z = [0.4, 0.6, 0.8, 1.0]
src_plane_num = len(src_z)
src_com_dist = [cosmos.comoving_distance(z).value * h for z in src_z]  # Mpc/h
for tag, z in enumerate(src_z):
    print("Source plane 1 at z = %.2f, %.5f Mpc/h" % (z, src_com_dist[tag]))

# show the shear/kappa map
bin_num = 100
Rmax = 1000
print("Max to %.3f Mpc/h" % (com_dist_len * Rmax / 3600 / 180 * numpy.pi))
if rank == 0:
    ra_bin = numpy.linspace(-Rmax, Rmax, bin_num + 1)
    dec_bin = numpy.linspace(-Rmax, Rmax, bin_num + 1)
    ra, dec = GGLensing_tool.get_pos(ra_bin, dec_bin, bin_num)

    shear_data = [GGLensing_tool.cal_shear(ra, dec, nfw, z) for z in src_z]

    img = Image_Plot(xpad=0.25, ypad=0.25)
    img.subplots(4, src_plane_num)

    for i in range(src_plane_num):
        kappa, g1, g2 = shear_data[i][0], shear_data[i][1], shear_data[i][2]
        g = numpy.sqrt(g1 ** 2 + g2 ** 2)

        img.scatter_pts(0, i, ra, dec, kappa)
        img.scatter_pts(1, i, ra, dec, g)
        img.scatter_pts(2, i, ra, dec, g1)
        img.scatter_pts(3, i, ra, dec, g2)

        img.axs[0][i].set_title("kappa from source plane at Z=%.3f" % src_z[i])
        img.axs[1][i].set_title("g from source plane at Z=%.3f" % src_z[i])
        img.axs[2][i].set_title("g1 from source plane at Z=%.3f" % src_z[i])
        img.axs[3][i].set_title("g2 from source plane at Z=%.3f" % src_z[i])

        for j in range(4):
            img.set_label(j, i, 0, "Dec. [arcsec]")
            img.set_label(j, i, 1, "R.A. [arcsec]")
    img.save_img("./pic/sample.png")
    img.close_img()
# img.show_img()
comm.Barrier()


# random position of source galaxies
total_num = 10000000
half_num = int(total_num / 2)
max_angle = Rmax  # arcsec
max_radius = com_dist_len * max_angle / 3600 * numpy.pi / 180
gal_x = numpy.random.uniform(-max_angle, max_angle, total_num).astype(dtype=numpy.float32)
gal_y = numpy.random.uniform(-max_angle, max_angle, total_num).astype(dtype=numpy.float32)
print(max_radius, " Mpc/h")

# for iz in range(src_plane_num):
shear_data = GGLensing_tool.cal_shear(gal_x, gal_y, nfw, src_z[rank]).astype(dtype=numpy.float32)
#     shear_datas.append(shear_data)

h5f = h5py.File("/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/shear_map_%d.hdf5" % rank, "w")
h5f["/ra"] = gal_x
h5f["/dec"] = gal_y
h5f["/z"] = numpy.array([src_z[rank]], dtype=numpy.float32)

h5f["/g1"] = shear_data[1]
h5f["/g2"] = shear_data[2]
h5f["/kappa"] = shear_data[0]
h5f["/theta"] = shear_data[3]

h5f.close()

plt_data = [[shear_data[0], numpy.sqrt(shear_data[1] ** 2 + shear_data[2] ** 2)],
            [shear_data[1], shear_data[2]]]

print("Shear/kappa from sources at z=%.2f " % src_z[rank])
img = Image_Plot(xpad=0.25, ypad=0.25)
img.subplots(2, 2)
for i in range(2):
    for j in range(2):
        img.scatter_pts(i, j, gal_x[:10000], gal_y[:10000], plt_data[i][j][:10000])
        img.set_label(i, j, 0, "Dec. [arcsec]")
        img.set_label(i, j, 1, "R.A. [arcsec]")
img.save_img("./pic/background_%d.png" % rank)
# img.show_img()
img.close_img()
