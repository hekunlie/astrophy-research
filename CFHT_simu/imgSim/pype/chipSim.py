## Generate simulated images with the mock catalogs constructed by mocksample.py

# Galaxy model: bulge + disk profile
# PSF: constant Gauss PSF with half light radius = 0.7 arcsec.
# Galaxy Shape: intrinsic + shear

# The parameters are detailed in configuration file 'voiceSim.param'

import galsim
import galsim.des
import imgMock as imk
import genStrFun as gsf
import pyfits as pf
import numpy as np
import time
from glob import glob
import logging
from collections import Counter
import os, sys

# the working directory
dirname, filename = os.path.split(os.path.abspath(__file__))
maindir = "/".join(dirname.split("/")[:-1]) + "/"

## read configuration file
print "^_^ Read configuration file ..."
confn = "defaultSim.param"
confnAbs = dirname + "/" + confn
params = gsf.read_param(confnAbs)

field = params["field"].lower()

if not field in ["cdfs1", "cdfs2", "cdfs3", "cdfs4"]:
    raise ValueError("!!! Parameter 'field' can only be cdfs1-4.")
seed_seed = int(field[-1])

##  directory
#maindir = os.path.dirname(os.getcwd()) + "/"
#maindir = params["dir_root"]
catdir  = maindir + params["dir_cat"]
figdir  = maindir + params["dir_fig"]
outdir  = maindir + params["dir_sim"]%field

## parameters
pixel_scale   = float(params["pixel_scale"])
nchip         = int(params["nchip"])
xmsize        = int(params["xsize_mosaic"])
ymsize        = int(params["ysize_mosaic"])
imrot         = float(params["pos_angle"]) * galsim.degrees
sigma_bkg     = float(params["bkg_sigma"])
zpoint        = float(params["zeropoint"])
epoch         = int(params["epoch"])
expxx         = int(params["expxx"])
chipxx        = int(params["chipxx"])
rotell        = float(params["rotate_ell"])
gain          = np.loadtxt(catdir + params["gain_file"])
expn          = catdir + params["dither_file"]%field
bkgsign       = catdir + params["bsig_file"]
shear_method  = params["shear_method"]
reduced_g1    = float(params["reduced_g1"])
reduced_g2    = float(params["reduced_g2"])
if "%s" in params["bsig_file"]: bkgsign = bkgsign%field

mcatn         = catdir + params["mock_catalog"]%field

#if expxx>5 : raise ValueError("!!! Parameter 'expxx' should be in [0,5]")
if chipxx>32: raise ValueError("!!! Parameter 'chipxx' should be in [0,32]")

# create the directory
irotell = int(rotell/45.0)
if irotell == 0:
    pass
elif irotell in [1, 2, 3]:
    outdir = outdir + "imgRot%d/"%int(rotell)
else:
    raise ValueError("!!! Please set a right 'rotate_ell' value.")

if os.path.exists(outdir):
    if len(glob(outdir+"*"))!=0: os.popen("rm -rf %s*"%outdir)
else:
    os.popen("mkdir %s"%outdir)

# log file
loggn = outdir + params["mock_catalog"].split(".")[0]%field + "img.log"
if os.path.exists(loggn): os.popen("rm %s"%loggn)
logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
logFile = logging.FileHandler(loggn)
logFile.setFormatter(logging.Formatter("%(message)s"))
logging.getLogger("chipsim_v3").addHandler(logFile)
logger = logging.getLogger("chipsim_v3")

if irotell in [1, 2, 3]:
    logger.info("")
    logger.info("    ** NOTE: You are rotating the intrinsic ellipticity by %.2f degrees."%rotell)
    logger.info("    **       It only works after the simulation with 'rotate_ell' = 0.0.")
    logger.info("")


logger.info("^_^ Simulating images for Field: %s."%field.upper())
## file name of each exposure and corresponding dither coordinates
exps = open(expn).read().splitlines()
nexp = len(exps)
# let me try to decide how many epochs by observed date in this field ...
pdate = []
for k in range(nexp):
    expk = exps[k]
    nkey = expk.split()[0]
    pdate += [nkey.split("T")[0]]

udate = Counter(pdate)
keys = udate.keys()
keys.sort()
epids = {}
nepoch = 0
for key in keys:
    nepoch += 1
    nk  = udate[key]
    id0 = pdate.index(key)
    epids[nepoch] = range(id0,id0+nk)

expids = range(nexp)
chipids = range(nchip)
logger.info("    Total %d epochs."%nepoch)
if epoch>=1:
    expids = epids[epoch]
    nepoch = 1
    logger.info("    Only simulating Epoch %d ..."%epoch)
    if expxx>=1:
        expids = [expids[expxx-1]]
        logger.info("    Only simulating Exposure %d ..."%(expxx))
if epoch==0 and expxx>=1:
    expids = []
    for kk in range(nepoch): expids += [epids[kk+1][expxx-1]]
    logger.info("    Only simulating Exposure %d ..."%(expxx))
if chipxx>=1:
    chipids = [chipids[chipxx-1]]
    logger.info("    Only simulating Chip %d ..."%chipxx)

logger.info("    Total %d epoch(s); %d exposure(s); %d chip(s) to be simulated."%(nepoch,len(expids),len(chipids)))

## background noise dispersion
if os.path.exists(bkgsign):
    logger.info("    Using external sigma file for background noise.")
    bkgsig = open(bkgsign).read().splitlines()[1:]
    nbkg = len(bkgsig)
    if nbkg != nexp: raise ValueError("!!! Exposures in %s is not equal to input."%params["bsig_file"]%field)
    bkgsig_dic = {}
    for i in range(nbkg):
        bkgsigxx = bkgsig[i].split()
        xn, xvalues = bkgsigxx[0], bkgsigxx[1:]
        bkgsig_dic[xn] = xvalues
else:
    logger.info("    The dispersion of background noise is fixed to %4.2f"%sigma_bkg)

## Gains
sid = np.argsort(gain[:,0])
gain = gain[sid]

## mock catalog
mcat  = pf.getdata(mcatn,ext=1)
nobj = mcat.size
ra, dec = mcat.field("ra"), mcat.field("dec")
mag = mcat.field("mag")
ell1, ell2 = mcat.field("e1"), mcat.field("e2")
ell, theta_e = mcat.field("e"), mcat.field("theta_e")
g1, g2 = mcat.field("g1"), mcat.field("g2")
g, theta_g = mcat.field("g"), mcat.field("theta_g")
hlr = mcat.field("scale_length")
bfrac = mcat.field("bulge_fraction")
ngal = sum(hlr>0.0)
npst = sum(hlr==0.0)
pid = np.arange(nobj)
pidt = pid[:ngal]
psid = pid[ngal:]

if shear_method=="constant":
    g1[:], g2[:] = reduced_g1, reduced_g2
    #g = np.sqrt(g1**2 + g2**2)
    gtt = galsim.Shear(g1=reduced_g1,g2=reduced_g2)
    g[:], theta_g[:] = gtt.g, gtt.beta.rad()*180.0/np.pi

logger.info("    Total %d objects in the mock catalog: %d galaxies and %d stars."%(nobj,ngal,npst))

## define the entire mosaic
logger.info("    Construct image mosaic of single exposure ...")
fimage = galsim.ImageF(xmsize, ymsize)
fimage.setOrigin(0,0)
imcen = fimage.trueCenter()
dudx = np.cos(imrot.rad()) * pixel_scale
dudy = -np.sin(imrot.rad()) * pixel_scale
dvdx = np.sin(imrot.rad()) * pixel_scale
dvdy = np.cos(imrot.rad()) * pixel_scale
affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=imcen)

tmpfile = "S2I3M4V5O6I7C8E"
title = "#ID PID X  Y  RA  DEC  mag bulge_frac  hlr e e1  e2  g g1  g2  beta_e beta_g  star\n"
fmt =   "%4d %7d %8.3f %8.3f %10.6f %10.6f %6.3f %5.3f %7.3f %5.3f %6.3f %6.3f %7.5f %8.5f %8.5f %7.3f %7.3f %2d\n"
logger.info("")
logger.info("^_^ Now begin the simulation ...")
t0 = time.time()
imglistn = outdir + "image.list"
imglist  = open(imglistn,"w")
iepoch = 0
for i in expids:
    iepoch += 1
    t1 = time.time()
    iexp = exps[i]
    ikey, ira, idec = iexp.split()
    logger.info("    Exp %s: (ra, dec)=(%s, %s)."%(ikey,ira,idec))
    sky_center = galsim.CelestialCoord(ra=float(ira)*galsim.degrees, dec=float(idec)*galsim.degrees)
    wcs = galsim.TanWCS(affine, sky_center, units=galsim.arcsec)
    fimage.wcs = wcs

    iexpn = outdir + ikey+".fits"
    #os.popen("cp %s %s"%(tmp_mos, iexpn))
    #os.popen("sethead %s RA=%s DEC=%s CRVAL1=%s CRVAL2=%s"%(iexpn,ira,idec,ira,idec))
    for j in chipids:
        logger.info("^_^ Simulating CHIP %d (EXP %d/%d)..."%(j+1,iepoch,len(expids)))
        jchipn = iexpn[:-5] + "_C%d.fits"%(j+1)
        x0, x1, y0, y1 = imk.chipDist_OmegaCam(j+1)
        jgain = gain[j,1]

        jbound = galsim.BoundsI(x0-1, x1-1, y0-1, y1-1)
        subimg = fimage.subImage(jbound)
        subimg.write(jchipn)
        comd1 = "sethead %s RA=%s DEC=%s EQUINOX=2000.0 RADECSYS=FK5 "
        comd2 = "IMAGEID=%d GAIN=%.5f ZP=%.2f SATLEVEL=%d AUTHOR=%s "
        comd3 = "OBJECT=VOICE_SIM"
        comd = comd1 + comd2 + comd3
        os.popen(comd%(jchipn,ira,idec,j+1,jgain,zpoint,65535,"DZL_PKU"))
        #print x0, x1, y0, y1

        # sky coverage of the single chip image
        os.popen("imsize -n 8 -dr %s > %s"%(jchipn,tmpfile))
        xxx = open(tmpfile).read().splitlines()
        os.popen("rm %s"%tmpfile)
        xx0 = xxx[0].split()
        xx1 = xxx[1].split()
        ra0, ra1 = float(xx0[2]), float(xx0[4])
        dec0, dec1 = float(xx1[1]), float(xx1[3])
        logger.info("    Sky coverage: RA %.4f - %.4f; DEC %.4f - %.4f."%(ra0,ra1,dec0,dec1))
        dra = (ra-ra0)*(ra-ra1)
        ddec = (dec-dec0)*(dec-dec1)
        pid0 = (dra<0.0) & (ddec<0.0) # corresponding to a mask array
        pi = list(pid[pid0])

        # Find image center
        hdr = pf.getheader(jchipn)
        satlevel = hdr["SATLEVEL"]
        xsize, ysize = hdr["NAXIS1"], hdr["NAXIS2"]
        os.popen("imsize -n 15 -d %s > %s"%(jchipn,tmpfile))
        xxx = open(tmpfile).read().splitlines()
        os.popen("rm %s"%tmpfile)
        xx0 = xxx[0].split()
        ra_cen, dec_cen = float(xx0[1]), float(xx0[2])

        logger.info("    Simulated image size %d*%d with pixel scale %.4f."%(xsize,ysize,pixel_scale))
        logger.info("    Image center: (Ra, Dec) = (%.4f %.4f)."%(ra_cen, dec_cen))
        logger.info("    Celestial projection type is TAN.")
        logger.info("    Total %d objects in the sky coverage."%len(pi))

        siminfon = jchipn[:-4]+"info"
        siminfo = open(siminfon,"w")
        siminfo.write(title)
        tbound = galsim.BoundsI(1,xsize,1,ysize)
        wcs = galsim.FitsWCS(jchipn)
        galimg = galsim.Image(wcs=wcs, bounds=tbound)

        kkk = 1
        psxxx, galxxx = 0, 0
        for k in pi:
            # first get xy positions through ra,dec
            rak, deck = ra[k], dec[k]
            world_pos = galsim.CelestialCoord(rak*galsim.degrees,deck*galsim.degrees)
            image_pos = galimg.wcs.toImage(world_pos)
            xk, yk = image_pos.x, image_pos.y

            if (xk-0.5)*(xk-xsize+0.5)>0.0 or (yk-0.5)*(yk-ysize+0.5)>0.0: continue

            # parameters for the object at image_pos
            magk   = mag[k]
            fluxk  = 10.0**(-0.4*(magk-zpoint))
            ek     = ell[k]
            theta_ek = theta_e[k]
            e1k    = ell1[k]
            e2k    = ell2[k]
            if magk<18.5 or magk>24.0: continue
            if abs(complex(e1k, e2k))>=1.0: continue

            if int(rotell/45.0)==1: e1k, e2k = -e2k,  e1k
            if int(rotell/45.0)==2: e1k, e2k = -e1k, -e2k
            if int(rotell/45.0)==3: e1k, e2k =  e2k, -e1k

            gk     = g[k]
            theta_gk = theta_g[k]
            g1k    = g1[k]
            g2k    = g2[k]
            hlrk   = hlr[k]*1.0
            bfrack = bfrac[k]
            #psf = galsim.Moffat(beta=6.0, fwhm=0.7, flux=1.0)
            psf = galsim.Gaussian(fwhm=0.7, flux=1.0)
            #psfzz = psf.drawImage()
            #psfzz.write("zzz.fits")
            #sys.exit(0)
            #if hlrk>0.0 and magk<18.5: continue

            # local coordinate and image offset
            x_nominal = xk + 0.5
            y_nominal = yk + 0.5
            ix_nominal = int(np.floor(x_nominal+0.5))
            iy_nominal = int(np.floor(y_nominal+0.5))
            dx = x_nominal - ix_nominal
            dy = y_nominal - iy_nominal
            offset = galsim.PositionD(dx,dy)
            local_wcs = galimg.wcs.local(image_pos)

            # model the object
            if k in psid:
                psk = 1
                # star model, i.e. PSF
                final = psf.withFlux(fluxk)
            else:
                psk = 0
                # galaxy model
                disk = galsim.Sersic(n=1.0,scale_radius=hlrk,flux=1.0)
                bulge = galsim.Sersic(n=4.0,half_light_radius=hlrk*0.6,flux=1.0)
                #if bfrack <= 0.5:
                #    disk = galsim.Sersic(n=1.0,scale_radius=hlrk,flux=1.0)
                #    bulge = galsim.Sersic(n=4.0,half_light_radius=hlrk*0.6,flux=1.0)
                #else:
                #    disk = galsim.Sersic(n=1.0,scale_radius=hlrk,flux=1.0)
                #    bulge = galsim.Sersic(n=4.0,half_light_radius=hlrk*0.6,flux=1.0)
                gal = bfrack*bulge+(1.0-bfrack)*disk
                gal = gal.withFlux(fluxk)
                gal_shape = galsim.Shear(g1=e1k, g2=e2k)
                #e1k, e2k, betak = gal_shape.e1, gal_shape.e2, gal_shape.beta.rad()*180.0/np.pi
                gal = gal.shear(gal_shape)
                gal_shape = galsim.Shear(g1=g1k, g2=g2k)
                #e1k, e2k, betak = gal_shape.e1, gal_shape.e2, gal_shape.beta.rad()*180.0/np.pi
                gal = gal.shear(gal_shape)
                final = galsim.Convolve([psf, gal])

            try:
                stamp = final.drawImage(wcs=local_wcs,offset=offset,method='no_pixel')
                stamp.setCenter(ix_nominal,iy_nominal)
                #if psk==0:
                #    if magk<=20.5:
                #        sxx = int(15.0*hlrk/pixel_scale)
                #        if np.mod(sxx,2)==1: sxx += 1
                #        sizeSub = max([sxx,10])
                #    elif magk<=24.5:
                #        sxx = int(12.0*hlrk/pixel_scale)
                #        if np.mod(sxx,2)==1: sxx += 1
                #        sizeSub = max([sxx,6])
                #    else: #elif magk<=25.0:
                #        sxx = int(12.0*hlrk/pixel_scale)
                #        if np.mod(sxx,2)==1: sxx += 1
                #        sizeSub = max([sxx,4])
                #    imageSub = galsim.ImageF(sizeSub, sizeSub)
                #    stamp = final.drawImage(image=imageSub,wcs=local_wcs,offset=offset,method='no_pixel')
                #    stamp.setCenter(ix_nominal,iy_nominal)

                #    fluxRatio = fluxk/stamp.added_flux
                #    stamp *= fluxRatio
                #else:
                #    stamp = final.drawImage(wcs=local_wcs,offset=offset,method='no_pixel')
                #    stamp.setCenter(ix_nominal,iy_nominal)

                bounds = stamp.bounds & galimg.bounds
                #if bounds.area() == 0.0: continue
                galimg[bounds] += stamp[bounds]

                sflux = np.sum(stamp.array)
                sizex, sizey = stamp.xmax-stamp.xmin+1, stamp.ymax-stamp.ymin+1
                #print
                oo = "    Flux ratio for Obj (type=%d) %4d: %5.3f; %4.2f mag; size (%3d,%3d)."
                sys.stdout.write("\r{0}".format(oo%(psk,kkk,fluxk/sflux,magk,sizex,sizey)))
                sys.stdout.flush()
                #logger.info(oo%(psk,kkk,fluxk/sflux,1./fluxRatio,magk,sizex,sizey))

                line = fmt%(kkk,k+1,xk,yk,rak,deck,magk,bfrack,hlrk,ek,e1k,e2k,gk,g1k,g2k,theta_ek,theta_gk,psk)
                siminfo.write(line)

                kkk += 1
                if psk==0: galxxx += 1
                if psk==1: psxxx  += 1
            except:
                pass

        logger.info("")
        siminfo.close()
        nbadgal = len(pi) - (psxxx + galxxx)
        logger.info("    Total %3d stars and %4d galaxies are projected; rejecting %d bad objects."%(psxxx,galxxx,nbadgal))

        # random Gaussian noise
        if os.path.exists(bkgsign):
            sigexp = bkgsig_dic[ikey]
            sigma_bkg = float(sigexp[j])
            logger.info("    Background Gaussian noise with dispersion of %5.2f."%sigma_bkg)
        if j+1<10: cbkgx = "0%s"%(j+1)
        if j+1>=10: cbkgx = "%s"%(j+1)
        seed_bkg = int(ikey[5:13]+ikey[14:]+cbkgx)
        rng = galsim.BaseDeviate(seed_bkg)
        gaussian_noise = galsim.GaussianNoise(rng, sigma=sigma_bkg)
        galimg.addNoise(gaussian_noise)
        galimg.write(jchipn)
        imgraw = pf.getdata(jchipn)
        imgraw[imgraw>65535.0]=65535.0
        if os.path.exists(jchipn): os.popen("rm %s"%jchipn)
        pf.writeto(jchipn,imgraw,hdr)
        imglist.write(jchipn+"\n")

    t2 = time.time()
    logger.info("^_^ Total time is %.2f minutes."%((t2-t1)/60.0))
    logger.info("")

imglist.close()
t3 = time.time()
logger.info("^_^ The simulation is done in %.2f minutes."%((t3-t0)/60.0))
os.popen("cp %s %s"%(confnAbs,outdir))
logger.info("^_^ The configuration file and simulated images are saved into %s."%params["dir_sim"]%field)

