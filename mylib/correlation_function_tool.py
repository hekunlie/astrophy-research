import tool_box
import numpy
from astropy.cosmology import FlatLambdaCDM
import camb
from camb import model, initialpower


def tomo_panel_num(tomo_bin_num):
    return int((tomo_bin_num * tomo_bin_num + tomo_bin_num) / 2)


def z2dist(zpts, H0, omega_m0):
    cosmos = FlatLambdaCDM(H0, omega_m0)
    return cosmos.comoving_distance(zpts).value * H0 / 100


def get_lrange(kmin, kmax, zmin, zmax, H0, Omeg_m0):
    cosmos = FlatLambdaCDM(H0=H0, Om0=Omeg_m0)
    dist_min = cosmos.comoving_distance(zmin).value * H0 / 100
    dist_max = cosmos.comoving_distance(zmax).value * H0 / 100
    lmax = kmax * dist_max
    lmin = kmin * dist_min
    print("%.4f < l <%.4f" % (lmin, lmax))
    return lmin, lmax


def get_nz(redshift, tomo_bin, redshift_e, redshift_e_bin_num, zlim=-1):
    # redshift_e is the expectation calculated from P(z)
    tomo_num = tomo_bin.shape[0] - 1
    # bin starts from 0
    if zlim > 0:
        ze_bin = numpy.linspace(0, zlim, redshift_e_bin_num + 1)
    else:
        ze_bin = numpy.linspace(0, redshift_e.max() + 0.0001, redshift_e_bin_num + 1)

    ze_bin_cent = (ze_bin[1:] + ze_bin[:-1]) / 2
    # one row for one tomo_bin
    zehist = numpy.zeros((tomo_num, redshift_e_bin_num))

    for i in range(tomo_num):
        idx1 = redshift >= tomo_bin[i]
        idx2 = redshift < tomo_bin[i + 1]
        idx = idx1 & idx2

        zz_e_num = numpy.histogram(redshift_e[idx], ze_bin)[0]
        zehist[i] = zz_e_num / zz_e_num.sum()

    return zehist, ze_bin, ze_bin_cent


def get_lenq(nz, com_dist):
    # nz: PDF of each tomo_bin, [tomo_bin_num, zhist_num],
    #       one row for one tomo_bin
    tomo_num, zpts_num = nz.shape
    qx = numpy.zeros((tomo_num, zpts_num))

    for i in range(tomo_num):
        nz_dist = nz[i] / com_dist
        for j in range(zpts_num):
            qx[i, j] = nz[i, j:].sum() - numpy.sum(nz_dist[j:]) * com_dist[j]
    return qx


def get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, zpts, kmax=3):
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


def get_PK(camb_result, kpts_num=300, kmin=1e-4, kmax=3):
    # camb_result: camb result instance
    kh_nonlin, z_nonlin, pk_nonlin = camb_result.get_matter_power_spectrum(minkh=kmin, maxkh=kmax, npoints=kpts_num)
    s8 = numpy.array(camb_result.get_sigma8())
    return z_nonlin, kh_nonlin, pk_nonlin, s8


def get_interpPk(cambPk, zpts, kpts):
    znum, knum = zpts.shape[0], kpts.shape[0]
    pk_result = numpy.zeros((znum, knum))
    for i in range(znum):
        pk_result[i] = cambPk.P(zpts[i], kpts)
    return pk_result


def get_PK_interp(camb_result):
    PK = camb_result.get_matter_power_interpolator(nonlinear=True, hubble_units=True, k_hunit=True)
    return PK


def ready4PL(Lmin, Lmax, Lpts_num, kmin, kmax, com_dist, zpts_num):
    # calculate the k corresponding to L (specified by user)
    Lpts = tool_box.set_bin_log(Lmin, Lmax, Lpts_num)
    integ_k = numpy.zeros((zpts_num, Lpts_num))

    dLpts = Lpts[1:] - Lpts[:-1]

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

    return Lpts, dLpts, integ_k


def get_PL(integ_pk, Lpts_num, integ_factor, delta_com_dist):
    PL = numpy.zeros((Lpts_num,))
    for i in range(Lpts_num):
        integ_part = integ_pk[:, i] * integ_factor
        PL[i] = numpy.sum((integ_part[1:] + integ_part[:-1]) / 2 * delta_com_dist)
    return PL