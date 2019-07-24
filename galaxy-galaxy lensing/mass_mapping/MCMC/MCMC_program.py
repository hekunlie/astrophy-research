import numpy


def gauss_weight(sig, radius):
    """weight for shear measurement"""
    return numpy.exp(-radius**2/2/sig**2)*40

def shear_profile(params, x, y):
    ''' the shear field '''
    ampl, dx, dy, sig = params[0], params[1], params[2], params[3]
    r = numpy.sqrt(x ** 2 + y ** 2)
    g = ampl * r/numpy.pi/2/sig/sig*numpy.exp(-((x - dx) ** 2 + (y - dy) ** 2) / 2/sig/sig)
    # if type(x) is not float:
    #     idx = numpy.abs(g) < 0.001
    #     g[idx] = 0
    #     idx = r == 0
    #     r[idx] = 1
    #     x[idx] = 0
    #     y[idx] = 0
    sin_theta = x/r
    cos_theta = y/r
    sin_2theta = 2*sin_theta*cos_theta
    cos_2theta = cos_theta**2 - sin_theta**2
    g1 = g*cos_2theta
    g2 = -g*sin_2theta
    return g, g1, g2

def shear_slope(params, x, y):
    ''' the shear field '''
    a1, a2, a3 = params[0], params[1], params[2]
    g = a1 + a2*x + a3*y
    return g

def ln_gh_prior(theta):
    tag = 0
    num = 4
    fit_range = [[0, 1.5],[-2, 2],[-2, 2],[1, 6]]
    for i in range(num):
        if fit_range[i][0] <= theta[i] <= fit_range[i][1]:
            tag += 1
    if tag == num:
        return 0.0
    else:
        return -numpy.inf

def ln_gh_prior_slope(theta):
    tag = 0
    num = 3
    fit_range = [[-0.1, 0.1],[-0.1, 0.1],[-0.1, 0.1]]
    for i in range(num):
        if fit_range[i][0] <= theta[i] <= fit_range[i][1]:
            tag += 1
    if tag == num:
        return numpy.exp(-(theta[0]**2+theta[1]**2+ theta[2]**2)/2/0.1/0.1)
    else:
        return -numpy.inf

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
        xi = numpy.sum(numpy.abs(g - fit_g)**2)*10000
        return lp - xi


def surface_fit(theta, pos, fxy):
    x, y = pos[0], pos[1]
    fun = shear_profile(theta, x, y)[0]
    return fun - fxy

def lst_g_fit(theta, pos, bins, bin_num2, inverse, tag):
    x, y, G, NU = pos[0], pos[1], pos[2], pos[3]
    g = shear_profile(theta, x, y)[tag]
    G_h = G - NU * g
    num = numpy.histogram(G_h, bins)[0]
    n1 = num[0:bin_num2][inverse]
    n2 = num[bin_num2:]
    xi = numpy.sum((n1 - n2) ** 2 / (n1 + n2)) * 0.5
    return xi


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
            # xi = numpy.sum(numpy.abs(n1 - n2))*0.5*10000
            return lp - xi
        return -numpy.inf

def ln_prob_g_slope(theta, G, NU, bins, bin_num2, inverse, x, y):
    """
    :param theta:
    :param G:
    :param NU:
    :param bins:
    :param bin_num2:
    :param inverse:
    :param x:
    :param y:
    :return:
    """
    lp = ln_gh_prior_slope(theta)
    if not numpy.isfinite(lp):
        return -numpy.inf
    else:
        g = shear_slope(theta, x, y)
        if numpy.abs(g).max() <= 0.1:
            G_h = G - NU * g
            num = numpy.histogram(G_h, bins)[0]
            n1 = num[0:bin_num2][inverse]
            n2 = num[bin_num2:]
            idx1 = n1 == 0
            idx2 = n2 == 0
            znum = numpy.sum(idx1 & idx2)
            if znum > 0:
                print("ZERO",znum)
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