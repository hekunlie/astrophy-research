def z2dist(zebin, H0, omega_m0, omega_k=0):
    cosmos = FlatLambdaCDM(H0, omega_m0)
    h = H0 / 100

    num = zebin.shape[0]
    dist_bin = numpy.zeros_like(zebin)

    for i in range(num):
        dist_bin[i] = cosmos.comoving_distance(zebin[i]).value * h

    return dist_bin


def get_nz(z, ztomo_bin, ze, ze_bin_num):
    ztomo_num = ztomo_bin.shape[0] - 1

    z_bin_num = int((z.max() - z.min()) / 0.01) + 2
    z_bin = [z.min() + 0.01 * i for i in range(z_bin_num)]

    ze_bin = numpy.linspace(ze.min(), ze.max(), ze_bin_num + 1)

    dist_bin = z2dist(ze_bin, 70, 0.27)

    zhist = []
    zehist = []

    img = Image_Plot()
    img.subplots(1, 4)

    for i in range(redshift_bin_num):
        idx1 = z >= ztomo_bin[i]
        idx2 = z < ztomo_bin[i + 1]
        idx = idx1 & idx2

        zz = z[idx]
        zz_e = ze[idx]

        zz_num = numpy.histogram(zz, z_bin)[0]
        zz_e_num = numpy.histogram(zz_e, ze_bin)[0]

        zz_e_dens = zz_e_num / zz_e_num.sum()

        zhist.append(zz_num)
        zehist.append(zz_e_dens)

        img.axs[0][0].hist(zz, z_bin, histtype="step", label="[%.2f, %.2f]" % (ztomo_bin[i], ztomo_bin[i + 1]))
        img.axs[0][1].plot(ze_bin[:-1], zz_e_dens, label="[%.2f, %.2f]" % (ztomo_bin[i], ztomo_bin[i + 1]))
        img.axs[0][2].plot(dist_bin[:-1], zz_e_dens, label="[%.3f, %.3f]" % (dist_bin[i], dist_bin[i + 1]))

    img.axs[0][3].plot(ze_bin, dist_bin)
    img.show_img()
    #     print(ze_bin)
    #     print(dist_bin)

    return zhist, z_bin, zehist, ze_bin


def get_lenq(dist, nz):
    nz_dist = nz / dist

    zpts_num = nz.shape[0]
    qx = numpy.zeros((zpts_num,))

    for i in range(zpts_num - 1):
        qx[i] = nz[i:].sum() - dist[i] * (nz_dist[i:].sum() - (nz_dist[i] + nz_dist[-1]) / 2)

    return qx


def get_pkl(H0, omega_m0, As, zmax, zbin, nz_i, nz_j):
    c = 299792.458  # Km/s
    h = H0 / 100  # Km/s/Mpc

    alpha = 9 / 4 * omega_m0 ** 2 * (H0 / c) ** 4  # Mpc^{-4}

    # midle point of each z bin
    zbin_mid = (zbin[1:] + zbin[:-1]) / 2
    zpts_num = zbin_mid.shape[0]

    scale_factor_sq_inve = (zbi_mid + 1) ** 2

    dist_bin_mid = z2dist(zbin_mid, H0, omega_m0)

    len_qi = get_lenq(dist_bin_mid, nz_i)
    len_qj = get_lenq(dist_bin_mid, nz_j)

    omg_b0h2 = 0.05 * (H0 / 100) ** 2
    omg_m0h2 = omega_m0 * (H0 / 100) ** 2

    # Now get matter power spectra and sigma8 at redshift 0 and 0.8
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=omg_b0h2, omch2=omg_m0h2)
    pars.InitPower.set_params(As=As, ns=0.965)
    # Note non-linear corrections couples to smaller scales than you want
    pars.set_matter_power(redshifts=zbin, kmax=3)

    # Non-Linear spectra (Halofit)
    pars.NonLinear = model.NonLinear_both
    results.calc_power_spectra(pars)
    kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=1e-4, maxkh=3, npoints=200)
    s8 = numpy.array(results.get_sigma8())