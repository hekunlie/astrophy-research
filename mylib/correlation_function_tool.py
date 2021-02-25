import tool_box
import numpy
from astropy.cosmology import FlatLambdaCDM
import camb
from camb import model, initialpower
import scipy



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


def get_nz(redshift, tomo_bin, redshift_e, redshift_e_bin_num, zlim=None, norm=True):
    # redshift_e is the expectation calculated from P(z)
    tomo_num = tomo_bin.shape[0] - 1
    # bin starts from 0
    if zlim:
        ze_bin = numpy.linspace(zlim[0], zlim[1], redshift_e_bin_num + 1)
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
        if norm:
            zehist[i] = zz_e_num / zz_e_num.sum()
        else:
            zehist[i] = zz_e_num

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
    sigma8 = results.get_sigma8()
    return results, sigma8, pars


def As2sigma8(As, Omega_cm0, Omega_bm0, zpts, H0, ns=0.965, interp_kmax=3):
    h = H0/100
    omg_cm0h2 = Omega_cm0*h*h
    omg_bm0h2 = Omega_bm0*h*h

    num = As.shape[0]
    sigma8 = numpy.zeros((num,))
    for i in range(num):
        s8 = get_CambResult(H0, omg_cm0h2[i], omg_bm0h2[i], As[i], ns, zpts, kmax=interp_kmax)[1]
        # print(s8,type(s8))
        sigma8[i] = s8[0]
    return sigma8

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

    # print("k h/Mpc ~[%.4f, %.4f]" % (integ_k.min(), integ_k.max()))
    numpy.savez("/mnt/perc/hklee/CFHT/correlation/test/integ_k.npz", integ_k)

    # only in the interval, [kmin, kmax], will be calculated
    idx = integ_k < kmin
    integ_k[idx] = -1
    idx = integ_k > kmax
    integ_k[idx] = -1
    numpy.savez("/mnt/perc/hklee/CFHT/correlation/test/integ_k_cut.npz", integ_k)
    return Lpts, dLpts, integ_k


def get_PL(integ_pk, Lpts_num, integ_factor, delta_com_dist):
    PL = numpy.zeros((Lpts_num,))
    for i in range(Lpts_num):
        integ_part = integ_pk[:, i] * integ_factor
        PL[i] = numpy.sum((integ_part[1:] + integ_part[:-1]) / 2 * delta_com_dist)
    return PL



def get_tomo_xi(As, Omega_cm0, Omega_bm0, h, zpts, inv_scale_factor_sq, zhist, z4pk_interp, theta_radian):
    '''
    calculate the xi_+/- using camb for mcmc or ...
    :param As: Amplitude of initial power spectrum
    :param Omega_cm0: cold dark matter density
    :param Omega_bm0: baryonic matter density
    :param h: H0 = 100*h Km/s/Mpc
    :param zpts: (n,) array, the redshift points for the integral, the centers of the redshift bins
    :param inv_scale_factor_sq: (n, ) array,  (1+zpts)^2
    :param zhist: (m, n) array, galaxy number density at each zpts point,
                    "m" tomographic z bin, "n" true z bin for each tomo z bin
    :param z4pk_interp: (n,), z points for power spectrum interpolation in camb,
                        [zmax, 0], at most 150 points
    :param theta_radian: (2m, n), theta, separation, in unit of radian
                        "m" tomographic z bin, "n" theta points
    :return: xi_+/- of each tomo panel at the given theta points,
            Pk(L) of each tomo panel at the given theta points,
            sigma8 at Z=0
    '''
    # t1 = time.time()
    ############### parameters  #############################
    # some of them should be modified for different purpose #
    #########################################################
    c = 299792.458  # Km/s
    ns = 0.965
    H0 = 100*h
    omg_cm0h2 = Omega_cm0 * h ** 2

    omg_bm0h2 = Omega_bm0 * h ** 2

    omg_m0 = Omega_cm0 + Omega_bm0

    alpha = 9 / 4 * omg_m0 ** 2 * (100 / c) ** 4  # h^4 \cdot Mpc^{-4}

    theta_num = theta_radian.shape[1]

    Lpts_min, Lpts_max, Lpts_num = 1, 15000, 8000

    kmin, cal_kmax, interp_kmax, kpts_num = 1e-4, 3, 4, 300

    # tomo panel num
    tomo_bin_num, zpts_num = zhist.shape
    tomo_panel_num = int((tomo_bin_num * tomo_bin_num + tomo_bin_num) / 2)

    # comoving distance Mpc/h
    com_dist = z2dist(zpts, H0, omg_m0)
    delta_com_dist = com_dist[1:] - com_dist[:-1]

    # lense efficiency
    qx = get_lenq(zhist, com_dist)

    # t2 = time.time()

    # g(x)^2/a(x)^2 in the integration
    integ_factor = numpy.zeros((tomo_panel_num, zpts_num))
    tag = 0
    for i in range(tomo_bin_num):
        for j in range(i, tomo_bin_num):
            integ_factor[tag] = qx[i] * qx[j] * inv_scale_factor_sq*alpha
            tag += 1
    # t3 = time.time()

    # set up bins for L of P(L), and decide the ks needed
    Lpts, dLpts, integ_kh = ready4PL(Lpts_min, Lpts_max, Lpts_num, kmin, cal_kmax, com_dist, zpts_num)
    numpy.savez("/mnt/perc/hklee/CFHT/correlation/test/com_dist.npz", com_dist,zpts)

    # L*Bessel_0(L*theta) in the final integral
    integ_Lpts_theta_p = numpy.zeros((tomo_panel_num, theta_num, Lpts_num))
    integ_Lpts_theta_m = numpy.zeros((tomo_panel_num, theta_num, Lpts_num))
    for i in range(tomo_panel_num):
        for j in range(theta_num):
            integ_Lpts_theta_p[i,j,:] = Lpts*scipy.special.j0(theta_radian[i,j] * Lpts)
            integ_Lpts_theta_m[i,j,:] = Lpts*scipy.special.jv(4,theta_radian[i,j] * Lpts)
    # t4 = time.time()

    # Pk interpolation
    integ_pk = numpy.zeros((zpts_num, Lpts_num))
    PLs = numpy.zeros((tomo_panel_num, Lpts_num))
    # the theoretical line
    xi_plus = numpy.zeros((tomo_panel_num, theta_num))
    xi_minus = numpy.zeros((tomo_panel_num, theta_num))
    xi_all = numpy.zeros((int(2*tomo_panel_num), theta_num))

    # calculate Pk using Camb
    camb_result = get_CambResult(H0, omg_cm0h2, omg_bm0h2, As, ns, z4pk_interp, kmax=interp_kmax)[0]
    pk_interp = camb_result.get_matter_power_interpolator(nonlinear=True, hubble_units=True, k_hunit=True)
    # sigma8 = camb_result.get_sigma8()
    # t5 = time.time()

    # get the Pk at given z points
    # only calculate the "k" < kmax ~ 2
    for i in range(zpts_num):
        idx = integ_kh[i] > 0
        integ_pk[i][idx] = pk_interp.P(zpts[i], integ_kh[i][idx])
    numpy.savez("/mnt/perc/hklee/CFHT/correlation/test/integ_pk.npz", integ_pk)
    # t6 = time.time()
    for i in range(tomo_panel_num):
        # get P(L)
        PLs[i] = get_PL(integ_pk, Lpts_num, integ_factor[i], delta_com_dist)

        # calculate the \chi_plus(\theta) , the theoretical line
        for j in range(theta_num):
            integ_part_p = integ_Lpts_theta_p[i,j,:] * PLs[i]
            xi_plus[i, j] = numpy.sum(((integ_part_p[1:] + integ_part_p[:-1]) / 2 * dLpts))

            integ_part_m = integ_Lpts_theta_m[i,j,:] * PLs[i]
            xi_minus[i, j] = numpy.sum(((integ_part_m[1:] + integ_part_m[:-1]) / 2 * dLpts))
    xi_all[:tomo_panel_num] = xi_plus
    # xi_all[tomo_panel_num:] = xi_minus
    # t7 = time.time()
    # print(t2-t1, t3-t2, t4-t3, t5-t4, t6-t5, t7-t6)
    return xi_plus/2/numpy.pi, xi_minus/2/numpy.pi, PLs, Lpts, xi_all/2/numpy.pi
    # return xi_plus/2/numpy.pi,PLs, Lpts# sigma8[-1],


