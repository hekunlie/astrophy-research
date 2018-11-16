# functions used for VOICE simulation

import os, sys
import numpy as np
import random
from astropy.io import fits
import galsim
from scipy.optimize import fmin
from scipy.interpolate import interp1d
import pylab as pl


def bfrac(n,seed=120000,figout=None):
    """
    Generate the distribution of bluge fraction
    """
    random.seed(seed)
    mu1, mu2 = 0.0, 1.0
    sig1, sig2 = 0.1, 0.001

    bfc = np.zeros(n)
    for i in range(n):
        p = random.uniform(0.0,1.0)
        if p < 0.9:
            ib = random.gauss(mu1,sig1)
            if ib<0.0: ib = abs(ib)
            if ib>1.0: ib = 1.0
        else:
            ib = random.gauss(mu2,sig2)
            if ib>1.0: ib = 1.0
            if ib<0.0: ib = 0.0

        bfc[i] = ib
    if figout is not None:
        pl.hist(bfc,bins=50,normed=True,histtype="step",color="black",linewidth=3)
        pl.xlabel("Bulge Fraction",fontsize=15)
        pl.savefig(figout)
        pl.clf()
        pl.close()

    return bfc


def e_disc(ne,seed=123000,figout=None):
    """
    Generate random ellipticity for a given number
    following the distribution of the ellipticity
    of disc-dominated galaxies

    See Miller et al, 2013, MNRAS
    """
    random.seed(seed)
    emin, emax = 0.0, 0.804
    a, e0 = 0.2539, 0.0256
    pe1 = lambda e: e*(1.0-np.exp((e-emax)/a))
    pe2 = lambda e: (1.0+e)*np.sqrt(e*e+e0*e0)
    pe  = lambda e: 2.4318*pe1(e)/pe2(e)

    pmin, pmax = 0.0, 2.02068
    # How to find pmin, pmax
    #numSteps = 10000000 # bigger the better but slower!
    #pmin = pe(xmin)
    #pmax = pmin
    #for i in range(numSteps):
    #    x = emin + (emax - emin) * float(i) / numSteps
    #    y = pe(x)
    #    if y < pmin: pmin = y
    #    if y > pmax: pmax = y
    #print pmin, pmax

    es = np.linspace(emin, emax, int((emax-emin)/0.000001)+1)
    pe_base = pe(es)
    pe_base = pe_base/np.sum(pe_base)

    rde = np.random.choice(es,ne,p=pe_base)

    #rde = np.zeros(ne)
    #es = np.linspace(emin, emax, int((emax-emin)/0.001)+1)
    #pe_base = pe(es)
    #for i in range(ne):
    #        x = random.uniform(emin,emax)
    #         y = random.uniform(pmin,pmax)
    #        while abs(y-pmax)<1.0e-3: y = random.uniform(pmin,pmax)
    #   if y > pe(x):
    #       xx = es[y<=pe_base]
    #            rde[i] = random.uniform(xx.min(), xx.max())
    #   else:
    #            rde[i] = x

    if figout is not None:
        eee = np.linspace(emin,emax,200)
        norm = np.trapz(pe(eee),eee)
        pl.plot(eee,pe(eee)/norm, color="red",linewidth=3,label="prior")
        pl.hist(rde,bins=50,normed=True,histtype="step",color="black",linewidth=3,label="Disc-dominated")
        pl.xlabel("ellipticity",fontsize=15)
        pl.legend()
        pl.savefig(figout)
        pl.clf()
        pl.close()
    return rde

def e_bulge(ne, seed=123400,figout=None):
    """
    Generate random ellipticity for a given number
    following the distribution of the ellipticity
    of bulge-dominated galaxies

    See Miller et al, 2013, MNRAS
    """
    random.seed(seed)
    b, c = 2.368, 6.691
    pe = lambda e: 27.7478*e*np.exp(-b*e - c*e*e)

    emin, emax = 0.0, 1.0
    pmin, pmax = 0.0, 2.64456

    es = np.linspace(emin, emax, int((emax-emin)/0.000001)+1)
    pe_base = pe(es)
    pe_base = pe_base/np.sum(pe_base)

    rbe = np.random.choice(es,ne,p=pe_base)

    #rbe = np.zeros(ne)
    #es = np.linspace(emin, emax, int((emax-emin)/0.001)+1)
    #pe_base = pe(es)
    #for i in range(ne):
    #        x = random.uniform(emin,emax)
    #        y = random.uniform(pmin,pmax)
    #        while abs(y-pmax)<1.0e-3: y = random.uniform(pmin,pmax)
    #        if y > pe(x):
    #       xx = es[y<=pe_base]
    #            rbe[i] = random.uniform(xx.min(), xx.max())
    #   else:
    #        rbe[i] = x

    if figout is not None:
        eee = np.linspace(emin,emax,200)
        norm = np.trapz(pe(eee),eee)
        pl.plot(eee,pe(eee)/norm, color="red", linewidth=3, label="prior")
        pl.hist(rbe,bins=50,normed=True,histtype="step",color="black",linewidth=3,label="Bulge-dominated")
        pl.xlabel("ellipticity",fontsize=15)
        pl.legend()
        pl.savefig(figout)
        pl.clf()
        pl.close()
    return rbe

def theta(n, seed=1234567):
    """
    Generate random orientation angle for galaxies
    """
    np.random.seed(seed)
    theta = np.random.uniform(-90.0, 90.0, n)

    return theta

def re_disk2(mag):
    """
    If one knows the magnitudes, then the scalelenghts can be derived directly.
    """
    rd = np.exp(-1.320-0.278*(mag-23.0)) # major-axis scalelength for r mag
    return rd

def re_disk_(mag, seed=123, figout=None):
    random.seed(seed)
    nobj = len(mag)
    rd = np.exp(-1.145 - 0.269 * (mag - 23.0))
    a, alpha = rd / 1.13, 4.0 / 3.0
    rr = np.zeros(nobj)
    rs = np.arange(0.035, 15.0, 0.001)

    abin = np.arange(a.min()-0.0005, a.max()+0.001, 0.001)
    nbin = len(abin) -1
    nxx = 0
    for i in range(nbin):
        a0, a1 = abin[i], abin[i+1]
        aMedian = 0.5*(a0+a1)
        idSub = (a > a0) & (a <= a1)
        nSub = idSub.sum()
        if nSub == 0: continue
        pre0 = rs*np.exp(-(rs/aMedian)**alpha)
        norm = np.sum(pre0)
        pre0 = pre0/norm
        rr[idSub] = np.random.choice(rs, nSub, p=pre0)
        nxx += nSub
    if nxx != nobj: raise ValueError("!!!nre and nmag are not identical.")

    if figout:
        pl.hist(rr, bins=80, normed=True, histtype="step", color="black", linewidth=3)
        pl.xlabel("Scalelength (arcsec)", fontsize=15)
        pl.xlim([0,3.5])
        pl.savefig(figout)
        pl.clf()
        pl.close()
    return rr


def re_disk(nre, mag, rlim=[0.1,5.0], seed=123450, figout=None):
    """
    Generate random scalelength for a given number
    following the distribution of the scalelength
    of dics-dominated galaxies. Unit: arcsec

    See Miller et al, 2013, MNRAS
    """
    random.seed(seed)
    nobj = len(mag)
    rmin, rmax = rlim[0], rlim[1]
    rs = np.linspace(rmin,rmax,int((rmax-rmin)/0.001)+1)
    rd = np.exp(-1.145-0.269*(mag-23.0)) # major-axis scalelength for i mag
    a, alpha = rd/1.13, 4.0/3.0
    spr = np.zeros(len(rs))
    # updated version
    for i in range(nobj):
        a0 = a[i]
        pre0 = rs*np.exp(-(rs/a0)**alpha)
        norm = np.trapz(pre0, rs)
        spr += pre0/norm
    spr = spr/float(nobj)
    pre = interp1d(rs, spr, kind="cubic")
    #spr = []
    #for i in range(nobj):
    #    a0 = a[i]
    #        pre0 = lambda r: r*np.exp(-(r/a0)**alpha)
    #        norm = np.trapz(pre0(rs), rs)
    #        pre1 = lambda r: pre0(r)/norm
    #
    #        spr += [pre1]
    #pre = lambda r: sum(i(r) for i in spr)/float(nobj)

    # How to find pmin, pmax
    #print "^_^ Looking for the max and min of the PDF ..."
    #r00 = fmin(lambda r: -pre(r), 0.2, disp=False)[0]
    #pmax = pre(r00)
    #pmin = min([pre(rmax), pre(rmin)])
    #print "   Max=%.4f; Min=%.4f"%(pmax, pmin)

    #print "^_^ Generating the random numbers ..."
    rs = np.linspace(rmin,rmax,int((rmax-rmin)/0.000001)+1)
    pre_base = pre(rs)
    pre_base = pre_base/np.sum(pre_base)

    rr = np.random.choice(rs,nre,p=pre_base)

    #rr = np.zeros(nre)
    #pre_base = pre(rs)
    #for j in range(nre):
    #        x = random.uniform(rmin,rmax)
    #   y = random.uniform(pmin,pmax)
    #        while abs(y-pmax)<1.0e-3: y = random.uniform(pmin,pmax)
    #        if y > pre(x):
    #       xx = rs[y<=pre_base]
    #            rr[j] = random.uniform(xx.min(), xx.max())
    #        else:
    #            rr[j] = x

    if figout is not None:
        rsxx = np.linspace(rmin,rmax,int((rmax-rmin)/0.01)+1)
        pl.plot(rsxx, pre(rsxx), color="red", linewidth=3, label="prior")
        pl.hist(rr,bins=150,normed=True,histtype="step",color="black",linewidth=3)
        pl.xscale("log")
        pl.xlabel("Scalelength (arcsec)",fontsize=15)
        pl.xlim([rmin, rmax])
        pl.legend()
        pl.savefig(figout)
        pl.clf()
        pl.close()
    #print "^_^ Done."
    return rr

def chipDist_OmegaCam(nchip):
    """
    calculate the edges in pixel for a given OmegaCam CCD chip
    """
    xt, yt = 17027, 16963
    x0, y0 = 2047, 4000
    g1 = 93
    g2 = 460
    g3 = 43

    mod = (nchip-1)/8
    xx0 = (mod+1)/2
    xx1 = mod - xx0

    echip = nchip - mod*8
    nx0 = 1 + (8-echip)*(x0+g1)
    nx1 = x0 + (8-echip)*(x0+g1)

    ny0 = 1 + (g2 + 4000)*xx0 + (g3 + 4000)*xx1
    ny1 = 4000 + (g2 + 4000)*xx0 + (g3 + 4000)*xx1

    return nx0, nx1, ny0, ny1

def eObs(e1,e2,g1,g2):
    """
    Calculate the sheared (observed) ellipticity using the
    intrinsic ellipticity and cosmic shear components.

    Parameters:
    e1, e2: scalar or numpy array
    g1, g2: scalar or numpy array

    Return:
    Sheared (observed) ellipticity components, absolute value, and orientation
    in format of scalar or numpy array

    ** NOTE: e1, e2, g1, and g2 should have the same dimensions.
    """
    if np.isscalar(e1):
        e1 = np.array([e1])
        e2 = np.array([e2])
        g1 = np.array([g1])
        g2 = np.array([g2])
    else:
        e1 = e1.flatten()
        e2 = e2.flatten()
        g1 = g1.flatten()
        g2 = g2.flatten()

    # calculate the sheared (observed) ellipticity using complex rule
    nobj = len(e1)
    e1obs, e2obs = [], []
    eeobs, theta = [], []
    for i in range(nobj):
        e = complex(e1[i], e2[i])
        g = complex(g1[i], g2[i])
        ee, gg = abs(e), abs(g)
        if gg<=1.0:
            tt = e + g
            bb = 1.0 + e*g.conjugate()
            eobs = tt/bb
        else:
            tt = 1.0 + g*e.conjugate()
            bb = e.conjugate() + g.conjugate()
            eobs = tt/bb

        # derive the orientation
        dd = 0.5*np.arctan(abs(eobs.imag/eobs.real))*180.0/np.pi
        if eobs.imag>0 and eobs.real>0: dd = dd
        if eobs.imag>0 and eobs.real<0: dd = 90.0 - dd
        if eobs.imag<0 and eobs.real>0: dd =  0.0 - dd
        if eobs.imag<0 and eobs.real<0: dd = dd - 90.0

        e1obs += [eobs.real]
        e2obs += [eobs.imag]
        eeobs += [abs(eobs)]
        theta += [dd]

    e1obs,e2obs,eeobs,theta = np.array(e1obs),np.array(e2obs),np.array(eeobs),np.array(theta)
    if nobj == 1: e1obs,e2obs,eeobs,theta = e1obs[0],e2obs[0],eeobs[0],theta[0]

    return e1obs, e2obs, eeobs, theta

def magBias(mag, bias):
    """
    Magnitude bias correction
    """
    mbfun = interp1d(mag, bias, kind="cubic")

    return mbfun

def fbin(data, nbins=4,lim=[10,100],bedge=None):
    """
    Given the bin number, find the number of data points and the data ID in each bin.

    Parameters:
    data: array, 1D
    nbins: number of bins

    Return:
    bins: list with length=nbins and contains the number of data points in each bin
    binid: dictionary
          the data ID in each bin
    """
    binid = {}
    if bedge is not None:
        edge = bedge
        nbins = len(bedge)-1
        for i in range(nbins):
            lim0, lim1 = bedge[i], bedge[i+1]
            id0 = np.where(np.logical_and(data>=lim0,data<lim1))[0]
            binid[i] = id0
    else:
        ndat = len(data)
        idnew = np.where(np.logical_and(data>=lim[0], data<=lim[1]))[0]
        nn = len(idnew)
        data0 = data[idnew]
        idd = np.argsort(data0)

        cc = nn/nbins
        bins = []
        for i in range(nbins): bins += [cc]
        if sum(bins)!=nn: bins[-1] = nn - cc*(nbins-1)

        edge = np.zeros(nbins+1)
        edge[0] = max([data0.min(), lim[0]])
        edge[-1] = min([data0.max(), lim[1]])

        nobj = 0
        for i in range(nbins):
            num = bins[i]
            if i<nbins-1:
                idd0 = idd[i*num:(i+1)*num]
                edge[i+1] = data0[idd[(i+1)*num]]
            else:
                idd0 = idd[-num:]
                edge[i+1] = data0[idd[-1]]
            binid[i] = idnew[idd0]
            nobj += len(idd0)
        print("^_^ Total %d objects; True %d"%(nobj, nn))

    return edge, binid

def bootstrap(sample, weight, nsamples=1000, conf=0.683):
    """
    Arguments:
    sample - input sample of values
    nsamples - number of samples to generate

    Returns bootstrap sample of statistic(computed by statfunc), bias and confidence interal.
    """
    # defin the statistical function
    statfun = lambda a, wgh: np.sum(a*wgh)/np.sum(wgh)
    msam = statfun(sample, weight)

    n = len(sample)
    X = []
    for i in range(nsamples):
        idd = np.random.choice(n,n)
        resample = sample[idd]
        reweight = weight[idd]
        x = statfun(resample, reweight)
        X.append(x)

    plower  = (1-conf)/2.0
    pupper  = 1 -plower
    symconf = (quantile(X, plower), quantile(X, pupper)) # corrected.
    errbar = 0.5*(symconf[1]-symconf[0])
    if (symconf[1]-msam)*(msam-symconf[0])<0:
        raise ValueError("!!! Something is wrong.")
    return errbar

def quantile(x, q, qtype=7):
    """
    Args:
    x - input data
    q - quantile
    qtype - algorithm

    Compute quantiles from input array x given q.For median,
    specify q=0.5.

    References:
    http://reference.wolfram.com/mathematica/ref/Quantile.html
    http://wiki.r-project.org/rwiki/doku.php?id=rdoc:stats:quantile

    Author:
    Ernesto P.Adorio Ph.D.
    UP Extension Program in Pampanga, Clark Field.
    """
    y = sorted(x)
    if not (1 <= qtype <= 9): return None  # error!

    # Parameters for the Hyndman and Fan algorithm
    abcd = [(0,   0, 1, 0), # inverse empirical distrib.function., R type 1
            (0.5, 0, 1, 0), # similar to type 1, averaged, R type 2
            (0.5, 0, 0, 0), # nearest order statistic,(SAS) R type 3
            (0,   0, 0, 1), # California linear interpolation, R type 4
            (0.5, 0, 0, 1), # hydrologists method, R type 5
            (0,   1, 0, 1), # mean-based estimate(Weibull method), (SPSS,Minitab), type 6
            (1,  -1, 0, 1), # mode-based method,(S, S-Plus), R type 7
            (1.0/3, 1.0/3, 0, 1), # median-unbiased ,  R type 8
            (3/8.0, 0.25, 0, 1)   # normal-unbiased, R type 9.
            ]

    a, b, c, d = abcd[qtype-1]
    n = len(x)
    g, j = np.modf( a + (n+b) * q -1)
    if j < 0:
        return y[0]
    elif j >= n:
        return y[n-1]   # oct. 8, 2010 y[n]???!! uncaught  off by 1 error!!!

    j = int(np.floor(j))
    if g ==  0:
        return y[j]
    else:
        return y[j] + (y[j+1]- y[j])* (c + d * g)

def imgExtract(ra, dec, inputImg="mosaic.fits", outimg="zz.fits"):
    """

    """
    raCen, decCen = 52.826383, -28.078758
    xmsize, ymsize = 37800, 37800
    xchip, ychip = 2047, 4000
    pixel_scale    = 0.213
    imrot          = 0.0 * galsim.degrees

    imcen = galsim.PositionD(x=xmsize/2.0+0.5, y=ymsize/2.0+0.5)
    dudx = np.cos(imrot.rad()) * pixel_scale
    dudy = -np.sin(imrot.rad()) * pixel_scale
    dvdx = np.sin(imrot.rad()) * pixel_scale
    dvdy = np.cos(imrot.rad()) * pixel_scale
    affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=imcen)

    sky_center = galsim.CelestialCoord(ra=float(raCen)*galsim.degrees, dec=float(decCen)*galsim.degrees)
    wcs = galsim.TanWCS(affine, sky_center, units=galsim.arcsec)
    #fimage.wcs = wcs

    world_pos = galsim.CelestialCoord(ra*galsim.degrees,dec*galsim.degrees)
    image_pos = wcs.toImage(world_pos)
    xk, yk = int(round(image_pos.x)), int(round(image_pos.y))

    x0, x1 = int(round(xk-1023.0)), int(round(xk+1023.0))
    y0, y1 = int(round(yk-1999.5)), int(round(yk+1999.5))

    #if os.path.exists(outimg): os.popen("rm %s"%outimg)
    comd = "imcopy %s[%d:%d,%d:%d] %s"
    os.popen(comd%(inputImg,x0,x1,y0,y1,outimg))

    imgcat = fits.getdata(outimg)

    if os.path.exists(outimg): os.popen("rm %s"%outimg)
    return imgcat

