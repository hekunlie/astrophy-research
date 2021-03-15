data_path = "/mnt/perc/hklee/Galaxy_Galaxy_lensing_test/cata/background/continue_source_z"
shear_tag = 0

delta_sigma = numpy.zeros((4, 24))

for i in range(24):
    h5f = h5py.File(data_path + "/sheared_para_%d.hdf5" % i, "r")
    ra = h5f["/ra"][()]
    dec = h5f["/dec"][()]
    z = h5f["/z"][()]
    com_dist = cosmos.comoving_distance(z).value * H_0 / 100

    crit_coeff = 1662895.2081868195 * com_dist / com_dist_len / (com_dist - com_dist_len) / (1 + len_z)

    src_pos = SkyCoord(ra=ra * units.arcsec, dec=dec * units.arcsec, frame="fk5")
    position_theta = len_pos.position_angle(src_pos).rad

    cos_theta = numpy.cos(position_theta)
    sin_theta = numpy.sin(position_theta)
    cos_2theta = numpy.cos(2 * position_theta)
    sin_2theta = numpy.sin(2 * position_theta)
    cos_4theta = numpy.cos(4 * position_theta)
    sin_4theta = numpy.sin(4 * position_theta)

    h5f = h5py.File(data_path + "/data/sheared_data_%d_noise_free.hdf5" % i, "r")
    data = h5f["/data"][()]
    h5f.close()

    mg1, mg2, mnr, mu, mv = data[:, 0], data[:, 1], data[:, 2], data[:, 3], data[:, 4]

    mg1r = (mg1 * cos_2theta - mg2 * sin_2theta) * crit_coeff
    mg2r = (mg1 * sin_2theta + mg2 * cos_2theta) * crit_coeff
    mur = mu * cos_4theta - mv * sin_4theta
    mnur1 = mnr + mur
    mnur2 = mnr - mur
    t1 = time.time()
    Ds, Ds_err = fq.find_shear(mg1r, mnur1, 10, left=-1000, right=1000, fit_num=20, chi_gap=100)[:2]
    Dsx, Dsx_err = fq.find_shear(mg1r, mnur1, 10, left=-1000, right=1000, fit_num=20, chi_gap=100)[:2]
    t2 = time.time()
    delta_sigma[:, i] = Ds, Ds_err, Dsx, Dsx_err
    print(delta_sigma[:, i], t2 - t1)

