from sys import path,argv
path.append("/home/hklee/work/mylib")
import hk_tool_box
import hk_gglensing_tool
import numpy
import h5py
import galsim
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord
from astropy import units
from mpi4py import MPI




# data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z"
data_path = "/home/hklee/work/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z_10"


cmd = int(argv[1])


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
Mass = 3*10 ** 13  # M_sun/h
conc = 6  # concentration
len_z = float(argv[2])#0.1  # redshift
halo_position = galsim.PositionD(0, 0)  # arcsec
com_dist_len = cosmos.comoving_distance(len_z).value * h  # Mpc/h
print("Lens plane at z = %.2f, %.5f Mpc/h" % (len_z, com_dist_len))

Rmin, Rmax = 0.05, 0.07 # Mpc/h
separation_min, separation_max = Rmin/com_dist_len/numpy.pi*180*3600, Rmax/com_dist_len/numpy.pi*180*3600


# lens profile
nfw = galsim.NFWHalo(Mass, conc, len_z, halo_position, omega_m0, omega_lam0)
CF = hk_gglensing_tool.Cosmos_flat(omega_m0, 100*h)
CF.NFW_profile_galsim((0,0), Mass, conc, len_z)

if cmd == 0:

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    numprocs = comm.Get_size()


    total_src_num = int(argv[3])  # 20000000
    fore_or_back = int(argv[4])

    # initialize the RandomState
    seed_ini = numpy.random.randint(1, 100000, 1)
    rng = numpy.random.RandomState(int(seed_ini + 1000 * rank))

    if rank == 0:
        print("Radius: %.3f Mpc/h ~ %.3f Mpc/h" % (Rmin, Rmax))
        print("Theta: %.3f arcsec ~ %.3f arcsec" % (separation_min, separation_max))

    # separation
    separation_ = numpy.abs(rng.normal(0, separation_max*1.5, total_src_num))
    separation_len = separation_max - separation_min
    idx = separation_ >= separation_len
    scale = numpy.divmod(separation_[idx], separation_len)[0]
    separation_[idx] = separation_[idx] - scale*separation_len
    separation = separation_ + separation_min

    print(separation.max(), separation.min())

    theta = rng.uniform(0, numpy.pi * 2, total_src_num)
    ra = separation * numpy.cos(theta)
    dec = separation * numpy.sin(theta)


    # magnitude & flux
    mag_s, mag_e = 21, 25.5
    mag = hk_tool_box.mag_generator(total_src_num, mag_s, mag_e).astype(dtype=numpy.float32)
    flux = hk_tool_box.mag_to_flux(mag).astype(dtype=numpy.float32)

    # galactic radius
    radius_s, radius_e = 0.35, 1.0

    radius = ((hk_tool_box.radii_from_mags(mag, radius_s, radius_e) + 0.8)/0.187).astype(dtype=numpy.float32)


    seed = rng.randint(1, 2000000000, int(total_src_num/10000))


    if fore_or_back == 0:
        # non-sheared source which have lower redshifts than the lens
        src_z = rng.uniform(0.05, len_z, total_src_num)

        h5f = h5py.File(data_path + "/params/non_sheared_para_%d.hdf5"%rank, "w")
        h5f["/z"] = src_z.astype(dtype=numpy.float32)
        h5f["/ra"] = ra.astype(dtype=numpy.float32)
        h5f["/dec"] = dec.astype(dtype=numpy.float32)

        h5f["/flux"] = flux.astype(dtype=numpy.float32)
        h5f["/radius"] = radius.astype(dtype=numpy.float32)

        h5f["/g1"] = numpy.zeros((total_src_num,), dtype=numpy.float32)
        h5f["/g2"] = numpy.zeros((total_src_num,), dtype=numpy.float32)
        h5f["/gamma1"] = numpy.zeros((total_src_num,), dtype=numpy.float32)
        h5f["/gamma2"] = numpy.zeros((total_src_num,), dtype=numpy.float32)

        h5f["/seed"] = seed
        h5f.close()

    elif fore_or_back == 1:
        # true background sources
        # src_z = len_z + 0.1 + numpy.abs(rng.normal(0, 0.25, total_src_num).astype(dtype=numpy.float32))
        # idx = src_z > 0.9
        # src_z[idx] = src_z[idx] - 0.6

        src_z = rng.uniform(len_z + 0.1, 0.9, total_src_num).astype(dtype=numpy.float32)

        shear_data = CF.get_shear(ra, dec, src_z).astype(dtype=numpy.float32)

        print("Z: ", src_z.max(), src_z.min())
        print("g: ", shear_data[1].min(), shear_data[1].max(), shear_data[2].min(), shear_data[2].max())
        h5f = h5py.File(data_path + "/params/sheared_para_%d.hdf5"%rank, "w")
        h5f["/z"] = src_z
        h5f["/ra"] = ra
        h5f["/dec"] = dec

        h5f["/flux"] = flux
        h5f["/radius"] = radius

        h5f["/kappa"] = shear_data[0]
        h5f["/gamma1"] = shear_data[1]
        h5f["/gamma2"] = shear_data[2]
        h5f["/g1"] = shear_data[3]
        h5f["/g2"] = shear_data[4]

        h5f["/seed"] = seed

        h5f.close()

    else:
        zmax = 2.0
        h5f = h5py.File("/home/hklee/work/DECALS/DECALS_shear_catalog_old/cat_inform/data_hist/hist.hdf5", "r")
        hist = h5f["/z_hist"][()][0]
        zbin = h5f["/z_bin"][()]
        h5f.close()
        zbin_mid = (zbin[1:] + zbin[:-1]) / 2
        idx = zbin_mid <= zmax
        zbin_mid = zbin_mid[idx]
        hist = hist[idx]

        zbin_new = numpy.linspace(0, zmax, 1001)
        zbin_new_mid = (zbin_new[1:] + zbin_new[:-1]) / 2

        hist_new = numpy.interp(zbin_new_mid, zbin_mid, hist)

        hist_new = hist_new / hist_new.sum()
        rn = numpy.linspace(0, zmax, 100000001)
        # src_z = numpy.random.choice(zbin_new_mid, total_src_num, p=hist_new)
        # src_z_m = numpy.random.choice(zbin_new_mid, total_src_num, p=hist_new)
        # src_z_m_err = numpy.random.normal(0, (1+src_z_m)*0.1)
        # src_z = src_z_m + src_z_m_err
        src_z = numpy.zeros((total_src_num, ))
        test_z = numpy.linspace(len_z + 0.05, zmax, 20)
        sub_num = int(total_src_num/20)
        for iz in range(len(test_z)):
            st, ed = int(iz*sub_num), int((iz+1)*sub_num)
            src_z[st:ed] = test_z[iz]

        shear_data = numpy.zeros((5, total_src_num), dtype=numpy.float32)

        idx = src_z > len_z
        # shear_data_ = hk_gglensing_tool.get_shear(nfw, ra[idx], dec[idx], src_z[idx])
        g1, g2 = CF.get_shear(ra[idx], dec[idx], src_z[idx])[:2]

        shear_data[1][idx] = g1
        shear_data[2][idx] = g2

        print("Z: ", src_z.max(), src_z.min())
        print("g: ", shear_data[1].min(), shear_data[1].max(), shear_data[2].min(), shear_data[2].max())
        h5f = h5py.File(data_path + "/params/sheared_para_%d.hdf5" % rank, "w")
        h5f["/z"] = src_z
        # h5f["/z_m"] = src_z_m
        h5f["/ra"] = ra
        h5f["/dec"] = dec

        h5f["/flux"] = flux
        h5f["/radius"] = radius

        # h5f["/kappa"] = shear_data[0]
        h5f["/gamma1"] = shear_data[1]
        h5f["/gamma2"] = shear_data[2]
        # h5f["/g1"] = shear_data[3]
        # h5f["/g2"] = shear_data[4]

        h5f["/seed"] = seed

        h5f.close()


else:
    separation_bin_num = int(argv[2])
    separation_bin = hk_tool_box.set_bin_log(Rmin, Rmax, separation_bin_num+1)
    print(len(separation_bin), separation_bin)

    src_nums = numpy.zeros((separation_bin_num,))

    data_type = ["noise_free", "noisy_cpp"]
    shear_type = ["sheared"]


    for st_tag, st in enumerate(shear_type):
        # read sheared parameters
        h5f = h5py.File(data_path + "/params/stack_%s_para.hdf5"%st, "r")
        src_z = h5f["/z"][()]
        src_ra = h5f["/ra"][()]
        src_dec = h5f["/dec"][()]
        src_g1 = h5f["/g1"][()]
        src_g2 = h5f["/g2"][()]
        h5f.close()

        # position and separation angle
        pos_len = SkyCoord(ra=0 * units.deg, dec=0 * units.deg, frame="fk5")
        pos_src = SkyCoord(ra=src_ra * units.deg, dec=src_dec * units.deg, frame="fk5")

        separation_radian = pos_len.separation(pos_src).radian
        separation_radius = separation_radian * com_dist_len

        print("Separation: ",separation_radius.min(), separation_radius.max(),src_ra.max())

        position_angle = pos_len.position_angle(pos_src).radian

        sin_2theta = numpy.sin(2 * position_angle)
        cos_2theta = numpy.cos(2 * position_angle)

        sin_4theta = numpy.sin(4 * position_angle)
        cos_4theta = numpy.cos(4 * position_angle)

        for dt_tag, dt in enumerate(data_type):
            h5f = h5py.File(data_path + "/data/stack_%s_%s.hdf5"%(st, dt),"r")
            src_data = h5f["/data"][()]
            h5f.close()

            mgt = src_data[:, 0] * cos_2theta - src_data[:, 1] * sin_2theta
            mgx = src_data[:, 0] * sin_2theta + src_data[:, 1] * cos_2theta

            mnu1 = src_data[:, 2] + src_data[:, 3] * cos_4theta - src_data[:, 4] * sin_4theta
            mnu2 = src_data[:, 2] - src_data[:, 3] * cos_4theta + src_data[:, 4] * sin_4theta

            h5f = h5py.File(data_path + "/data/segment_%s_%s.hdf5"%(st, dt),"w")

            for ir in range(separation_bin_num):
                idx_r1 = separation_radius >= separation_bin[ir]
                idx_r2 = separation_radius < separation_bin[ir+1]
                idx = idx_r1 & idx_r2

                h5f['/%d/mgt'%ir] = mgt[idx]
                h5f['/%d/mgx'%ir] = mgx[idx]
                h5f['/%d/mnu1'%ir] = mnu1[idx]
                h5f['/%d/mnu2'%ir] = mnu2[idx]

                h5f['/%d/z'%ir] = src_z[idx]
                h5f['/%d/ra'%ir] = src_ra[idx]
                h5f['/%d/dec'%ir] = src_dec[idx]
                h5f['/%d/g1'%ir] = src_g1[idx]
                h5f['/%d/g2'%ir] = src_g2[idx]
                h5f['/%d/sep_radius'%ir] = separation_radius[idx]
                h5f['/%d/sep_radian'%ir] = separation_radian[idx]

                src_nums[ir] = idx.sum()

            h5f["/src_num"] = src_nums.astype(numpy.intc)
            h5f.close()

            print(st, dt, src_nums)

