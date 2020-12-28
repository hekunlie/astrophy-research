def get_factorial(x):
    if x == 0:
        return 1.
    else:
        j = 1.
        for i in range(1, x + 1):
            j *= i
        return j


def bessel_zero(x, m):
    bs = 1.
    one = 1
    x2 = x / 2.

    for i in range(1, m):
        one = -one
        a = x2 ** (i)
        b = 1.
        for j in range(1, i + 1):
            b *= j
        bs += one * (a / b) ** 2
    return bs


def z2dist(zpts, H0, omega_m0, omega_k=0):
    cosmos = FlatLambdaCDM(H0, omega_m0)
    h = H0 / 100

    com_dist = cosmos.comoving_distance(zpts).value * h
    delta_com_dist = com_dist[1:] - com_dist[:-1]

    return com_dist, delta_com_dist


def get_lrange(kmin, kmax, zmin, zmax, H0, Omeg_m0):
    cosmos = FlatLambdaCDM(H0=H0, Om0=Omeg_m0)
    dist_min = cosmos.comoving_distance(zmax).value * H0 / 100
    dist_max = cosmos.comoving_distance(zmax).value * H0 / 100
    lmax = kmax * dist_max
    lmin = kmin * dist_min
    print("%.4f < l <%.4f" % (lmin, lmax))
    return lmin, lmax


def get_nz(z, ztomo_bin, ze, ze_bin_num):
    ztomo_num = ztomo_bin.shape[0] - 1

    z_bin_num = int((z.max() - z.min()) / 0.01) + 2
    z_bin = [z.min() + 0.01 * i for i in range(z_bin_num)]

    ze_bin = numpy.linspace(0, ze.max() + 0.0001, ze_bin_num + 1)
    #     # add a zero to the head of ze_bin for the integral later starts from 0.
    #     ze_bin_cent = numpy.zeros_like(ze_bin)
    ze_bin_cent = (ze_bin[1:] + ze_bin[:-1]) / 2

    zhist = []
    zehist = []

    img = Image_Plot(xpad=0.25)
    img.subplots(1, 2)

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
        img.axs[0][1].plot(ze_bin_cent, zz_e_dens, label="[%.2f, %.2f]" % (ztomo_bin[i], ztomo_bin[i + 1]))

    img.show_img()

    return zhist, z_bin, zehist, ze_bin, ze_bin_cent


def get_lenq(com_dist, nz):
    zpts_num = nz.shape[0]
    qx = numpy.zeros((zpts_num,))
    nz_dist = nz / com_dist

    for i in range(zpts_num):
        qx[i] = nz[i:].sum() + numpy.sum(nz_dist[i:]) * com_dist[i]

    return qx


def get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, zpts, kpts_num=300, kmin=1e-4, kmax=2):
    # set parameters
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=omg_bm0h2, omch2=omg_cm0h2)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_matter_power(redshifts=zpts, kmax=kmax)

    # Non-Linear spectra (Halofit)
    pars.NonLinear = model.NonLinear_both
    pars.NonLinearModel.set_params(halofit_version='takahashi')
    # get results, it must run before the power speactra calculation
    results = camb.get_results(pars)

    return results, pars


def get_PK(H0, omg_cm0h2, omg_bm0h2, As, ns, zpts, kpts_num=300, kmin=1e-4, kmax=3):
    # set parameters
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=omg_bm0h2, omch2=omg_cm0h2)
    pars.InitPower.set_params(As=As, ns=ns)
    pars.set_matter_power(redshifts=zpts, kmax=kmax)

    # Non-Linear spectra (Halofit)
    pars.NonLinear = model.NonLinear_both
    pars.NonLinearModel.set_params(halofit_version='takahashi')
    # get results, it must run before the power speactra calculation
    results = camb.get_results(pars)

    results.calc_power_spectra(pars)
    kh_nonlin, z_nonlin, pk_nonlin = results.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints=kpts_num)
    s8 = numpy.array(results.get_sigma8())

    return z_nonlin, kh_nonlin, pk_nonlin, s8


def get_interpPk(cambPk, zpts, kpts):
    znum, knum = zpts.shape[0], kpts.shape[0]
    pk_result = numpy.zeros((znum, knum))
    for i in range(znum):
        pk_result[i] = cambPk.P(zpts[i], kpts)
    return pk_result


def get_PK_interp(H0, omg_cm0h2, omg_bm0h2, As, ns, zpts, kmax=5):
    # set parameters
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=H0, ombh2=omg_bm0h2, omch2=omg_cm0h2)
    pars.InitPower.set_params(As=As, ns=ns)
    # Non-Linear spectra (Halofit)
    pars.NonLinearModel.set_params(halofit_version='takahashi')

    PK = camb.get_matter_power_interpolator(pars, nonlinear=True, hubble_units=True, k_hunit=True, kmax=kmax, zs=zpts)

    return PK


def ready4PL(Lmin, Lmax, Lpts_num, kmin, kmax, com_dist, zpts_num):
    # calculate the k corresponding to L (specified by user)
    Lpts = tool_box.set_bin_log(Lmin, Lmax, Lpts_num)
    integ_k = numpy.zeros((zpts_num, Lpts_num))

    for i in range(Lpts_num):
        # integrate PK for each fixed L
        # k = L/com_dist at Z
        integ_k[:, i] = Lpts[i] / com_dist

    print("k h/Mpc ~[%.4f, %.4f]" % (integ_k.min(), integ_k.max()))

    # only in the interval, [kmin, kmax], will be calculated
    idx = integ_k < kmin
    integ_k[idx] = 0
    idx = integ_k > kmax
    integ_k[idx] = 0

    return Lpts, integ_k


def get_PL(pk_integ, factor_integ, delta_com_dist):
    Lpts_num = pk_integ.shape[1]
    PL = numpy.zeros((Lpts_num,))

    for i in range(Lpts_num):
        integ_part = pk_integ[:, i] * factor_integ
        PL[i] = numpy.sum((integ_part[1:] + integ_part[:-1]) / 2 * delta_com_dist)

    return PL