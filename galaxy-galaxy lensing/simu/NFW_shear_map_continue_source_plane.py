from sys import path,argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import gglensing_tool
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

data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z"

total_src_num = int(argv[1])#20000000
fore_or_back = int(argv[2])

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

Rmin, Rmax = 0.04, 10 # Mpc/h
separation_min, separation_max = Rmin/com_dist_len/numpy.pi*180*3600, Rmax/com_dist_len/numpy.pi*180*3600
separation_bin_num = numprocs
separation_bin = numpy.linspace(separation_min, separation_max, numprocs+1)
if rank == 0:
    print("Radius: %.3f Mpc/h ~ %.3f Mpc/h"%(Rmin, Rmax))
    print("Theta: %.3f arcsec ~ %.3f arcsec"%(separation_min, separation_max))

# initialize the RandomState
seed_ini = numpy.random.randint(1, 100000, 1)
rng = numpy.random.RandomState(int(seed_ini + 1000 * rank))


# lens profile
nfw = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)

# separation
separation = rng.uniform(separation_bin[rank], separation_bin[rank + 1], total_src_num)
theta = rng.uniform(0, numpy.pi * 2, total_src_num)
ra = separation * numpy.cos(theta)
dec = separation * numpy.sin(theta)

seed = rng.randint(1, 2000000000, int(total_src_num/10000))

if fore_or_back == 0:
    # non-sheared source which have lower redshifts than the lens
    src_z = rng.uniform(0.05, len_z, total_src_num)

    h5f = h5py.File(data_path + "/non_sheared_para_%d.hdf5"%rank, "w")
    h5f["/z"] = src_z.astype(dtype=numpy.float32)
    h5f["/ra"] = ra.astype(dtype=numpy.float32)
    h5f["/dec"] = dec.astype(dtype=numpy.float32)
    h5f["/g1"] = numpy.zeros((total_src_num,), dtype=numpy.float32)
    h5f["/g2"] = numpy.zeros((total_src_num,), dtype=numpy.float32)
    h5f["/kappa"] = numpy.zeros((total_src_num,), dtype=numpy.float32)
    h5f["/seed"] = seed
    h5f.close()

else:
    # true background sources
    src_z = len_z + 0.1 + numpy.abs(rng.normal(0, 0.35, total_src_num).astype(dtype=numpy.float32))

    shear_data = gglensing_tool.cal_shear(ra, dec, nfw, src_z).astype(dtype=numpy.float32)

    h5f = h5py.File(data_path + "/sheared_para_%d.hdf5"%rank, "w")
    h5f["/z"] = src_z
    h5f["/ra"] = ra
    h5f["/dec"] = dec

    h5f["/g1"] = shear_data[1]
    h5f["/g2"] = shear_data[2]
    h5f["/kappa"] = shear_data[0]
    h5f["/theta"] = shear_data[3]

    h5f["/seed"] = seed
    h5f.close()

