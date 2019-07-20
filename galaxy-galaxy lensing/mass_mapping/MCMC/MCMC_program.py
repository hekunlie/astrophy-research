import numpy



def ln_gh_prior(theta):
    tag = 0
    num = 4
    fit_range = [[0, 2],[-10, 10],[-10, 10],[1, 25]]
    for i in range(num):
        if fit_range[i][0] <= theta[i] <= fit_range[i][1]:
            tag += 1
    if tag == num:
        return 0.0
    else:
        return -numpy.inf


def shear_profile(params, x, y):
    ''' the shear field '''
    ampl, dx, dy, sig = params[0], params[1], params[2], params[3]
    r = numpy.sqrt(x ** 2 + y ** 2)
    g = ampl * r/numpy.pi/2/sig/sig*numpy.exp(-((x - dx) ** 2 + (y - dy) ** 2) / 2/sig/sig)
    if type(x) is not float:
        idx = numpy.abs(g) < 0.001
        g[idx] = 0
        idx = r == 0
        r[idx] = 1
        x[idx] = 0
        y[idx] = 0
    sin_theta = x/r
    cos_theta = y/r
    sin_2theta = 2*sin_theta*cos_theta
    cos_2theta = cos_theta**2 - sin_theta**2
    g1 = g*cos_2theta
    g2 = -g*sin_2theta
    return g, g1, g2


def ln_prob_g_test(theta, g, x, y, tag):
    """
    :param theta:
    :param g:
    :param x:
    :param y:
    :param tag: 0: g, 1: g1, 2: g2
    :return:
    """
    lp = ln_gh_prior(theta)
    if not numpy.isfinite(lp):
        return -numpy.inf
    else:
        fit_g = shear_profile(theta, x, y)[tag]
        xi = numpy.sum(numpy.abs(g - fit_g))
        return lp - xi


def surface_fit(theta, pos, fxy, tag):
    x, y = pos[0], pos[1]
    fun = shear_profile(theta, x, y)[tag]
    return fun - fxy


def ln_prob_g(theta, G, NU, bins, bin_num2, inverse, x, y, tag):
    """
    tag = 1: for g1, G=G1, NU = N + U
    tag = 2: for g2, G=G2, NU = N - U
    :param theta:
    :param G:
    :param NU:
    :param bins:
    :param bin_num2:
    :param inverse:
    :param x:
    :param y:
    :param tag:
    :return:
    """
    lp = ln_gh_prior(theta)
    if not numpy.isfinite(lp):
        return -numpy.inf
    else:
        g = shear_profile(theta, x, y)[tag]
        if numpy.abs(g).max() <= 0.1:
            G_h = G - NU * g
            num = numpy.histogram(G_h, bins)[0]
            n1 = num[0:bin_num2][inverse]
            n2 = num[bin_num2:]
            xi = numpy.sum((n1 - n2)**2/(n1 + n2))*0.5
            return lp - xi
        return -numpy.inf



def result_fun(params, ra, dec, radius):
    f_sub = radius/ numpy.pi / 2 / params[3][0] * numpy.exp(-((dec - params[1][0]) ** 2 + (ra - params[2][0]) ** 2) / 2 / params[3][0])
    f = params[0][0] * f_sub
    f_sig = f_sub*(params[0][2]+\
            params[0][0]*(((ra-params[1][0])**2 +(dec-params[2][0])**2)/2/params[3][0]**2 - 1/params[3][0])*params[3][3]+\
            params[0][0]/params[3][0]*(ra-params[1][0])*params[1][3]+\
            params[0][0]/params[3][0]*(dec-params[2][0])*params[2][3])
    return f, f_sig