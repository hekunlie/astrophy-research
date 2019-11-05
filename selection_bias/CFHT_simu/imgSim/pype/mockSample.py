# Generate mock catalogs for VOICE simulation

# Catalogs are generated tile by tile. There are four tiles for VOICE survey,
# so four mock catalogs will be generated.

# The parameters are detailed in configuration file 'voiceSim.param'

import galsim
from astropy.io import fits
import numpy as np
import time
import imgMock as imk
import genStrFun as gsf
import logging
import os, sys

t0 = time.time()

# the working directory
dirname, filename = os.path.split(os.path.abspath(__file__))
maindir = "/".join(dirname.split("/")[:-1]) + "/"

## read configuration file
print "^_^ Read configuration file."
confn = "defaultSim.param"
confnAbs = dirname + "/" + confn
params = gsf.read_param(confnAbs)

field = params["field"].lower()

if not field in ["cdfs1", "cdfs2", "cdfs3", "cdfs4"]:
    raise ValueError("!!! Parameter 'field' can only be 'cdfs1', 'cdfs2', 'cdfs3' or 'cdfs4'")
seed_seed = int(field[-1])

##  directory
#maindir = os.path.dirname(os.getcwd()) + "/"
#maindir = params["dir_root"]
catdir  = maindir + params["dir_cat"]
figdir  = maindir + params["dir_fig"]

## parameters
pixel_scale   = float(params["pixel_scale"])
pos_th     = float(params["pos_angle"]) * galsim.degrees
xmsize        = int(params["xsize_mosaic"])
ymsize        = int(params["ysize_mosaic"])
shear_method  = params["shear_method"]
reduced_g1    = float(params["reduced_g1"])
reduced_g2    = float(params["reduced_g2"])
ra_cen_ps     = float(params["ra_ps"])
dec_cen_ps    = float(params["dec_ps"])
ps_seed       = int(params["seed_ps"])
ps_file       = catdir + params["ps_file"]
grid_spacing  = float(params["grid_spacing"])
shear_catn    = catdir + params["shear_cat"]
if "%s" in params["shear_cat"]: shear_catn = shear_catn%field

galn          = catdir + params["galaxy_refcat"]%field
pstn          = catdir + params["star_refcat"]%field
mockn         = catdir + params["mock_catalog"]%field
contn         = catdir + params["catalog_contents"]
keyn          = params["mock_catalog"].split(".")[0]%field

seed_bfrac    = int(params["seed_bulge_fraction"]) + 10*2**seed_seed
seed_edisc    = int(params["seed_ell_disc"]) + 10*2**seed_seed
seed_ebulge   = int(params["seed_ell_bulge"]) + 10*2**seed_seed
seed_hdisc    = int(params["seed_hlr_disc"]) + 10*2**seed_seed
seed_hbulge   = int(params["seed_hlr_bulge"]) + 10*2**seed_seed
seed_theta    = int(params["seed_theta"]) + 10*2**seed_seed

if not shear_method in ["constant", "ps_shear", "extra"]:
    raise ValueError("!!! Please set a right 'shear_method' parameter.")

# log file
loggn = catdir + keyn + "sample.log"
if os.path.exists(loggn): os.popen("rm %s"%loggn)
logging.basicConfig(format="%(message)s", level=logging.INFO, stream=sys.stdout)
logFile = logging.FileHandler(loggn)
logFile.setFormatter(logging.Formatter("%(message)s"))
logging.getLogger("mocksample").addHandler(logFile)
logger = logging.getLogger("mocksample")

# field to be simulated
logger.info( "^_^ Generate the mock catalog for Field: %s."%field.upper())
logger.info("^_^ Now begin the simulation ...")
galcat = fits.getdata(galn,ext=1)
ra_gal, dec_gal = galcat.field("ra"), galcat.field("dec")
mag_gal   = galcat.field("mag")
ngal = galcat.size

pstcat = fits.getdata(pstn,ext=1)
ra_pst, dec_pst = pstcat.field("ra"), pstcat.field("dec")
mag_pst   = pstcat.field("mag")
npst = pstcat.size

nobj = ngal + npst
ra, dec = np.append(ra_gal,ra_pst), np.append(dec_gal,dec_pst)
mag = np.append(mag_gal,mag_pst)
logger.info("    Total %d galaxies and %d stars."%(ngal,npst))

# galaxy parameters
index = np.arange(1,nobj+1)
ell, theta, hlr, bfrac = np.zeros(nobj), np.zeros(nobj), np.zeros(nobj), np.zeros(nobj)
ell1, ell2 = np.zeros(nobj), np.zeros(nobj)
g1, g2 = np.zeros(nobj), np.zeros(nobj)
gt, beta = np.zeros(nobj), np.zeros(nobj)
eg1, eg2 = np.zeros(nobj), np.zeros(nobj)
egt, egd = np.zeros(nobj), np.zeros(nobj)

## galaxy background shear construction
if shear_method == "ps_shear":
    logger.info("    Background shear field construction ...")
    ps_image_xsize = int(xmsize*np.cos(pos_th.rad())+ymsize*np.sin(pos_th.rad()))
    ps_image_ysize = int(xmsize*np.sin(pos_th.rad())+ymsize*np.cos(pos_th.rad()))
    ps_image_size  = max([ps_image_xsize, ps_image_ysize])
    ps_image_size  = int(2.5*ps_image_size)
    ps_image = galsim.ImageF(ps_image_size, ps_image_size)
    ps_image.setOrigin(0,0)
    dudx = np.cos(pos_th.rad()) * pixel_scale
    dudy = -np.sin(pos_th.rad()) * pixel_scale
    dvdx = np.sin(pos_th.rad()) * pixel_scale
    dvdy = np.cos(pos_th.rad()) * pixel_scale
    affine = galsim.AffineTransform(dudx, dudy, dvdx, dvdy, origin=ps_image.trueCenter())
    sky_center = galsim.CelestialCoord(ra=ra_cen_ps*galsim.degrees, dec=dec_cen_ps*galsim.degrees)
    wcs = galsim.TanWCS(affine, sky_center, units=galsim.arcsec)
    ps_image.wcs = wcs

    rng_ps = galsim.BaseDeviate(ps_seed)
    ps = galsim.PowerSpectrum(ps_file, units = galsim.radians)
    ps.buildGrid(grid_spacing = grid_spacing,
                 ngrid = int(1+np.ceil(ps_image_size*pixel_scale/grid_spacing)),
                 rng = rng_ps.duplicate())
elif shear_method == 'constant':
    reduced_shear = np.sqrt(reduced_g1**2 + reduced_g2**2)
    logger.info("    Background shear field is assumed to be constant |g|=%.4f"%reduced_shear)
    g1[:ngal], g2[:ngal] = reduced_g1, reduced_g2
    #gt[:ngal] = reduced_shear
    #beta[:ngal] = imk.theta(ngal,seed=beta_seed)
else: # shear_method=="extra":
    logger.info("    Background shear is from extra catalog %s."%shear_catn.split("/")[-1])
    if not os.path.exists(shear_catn):
        raise ValueError("!!! Put the shear catalog to %s."%catdir)
    shearCat = np.loadtxt(shear_catn)
    nshear = shearCat.shape[0]
    if nshear != ngal:
        raise ValueError("!!! Shear number (%d) is NOT equal to input galaxy number (%d)."%(nshear,ngal))
    g1[:ngal], g2[:ngal] = shearCat[:,0], shearCat[:,1]

logger.info("")
logger.info("    Generate galaxy shape and profile parameters ... ")
logger.info("    Galaxy model priors: bulge-fraction, ellipticity, half-light radius ...")
#bfracGal = imk.bfrac(ngal,seed=seed_bfrac,figout=figdir+keyn+"_bulgeFrac.png")
bfracGal = imk.bfrac(ngal,seed=seed_bfrac,figout=None)
id0, id1 = np.where(bfracGal<=0.5)[0], np.where(bfracGal>0.5)[0]
ndisc = len(id0)
nbulge = len(id1)
logger.info("    Total %d disc-dominated and %d bulge-dominated galaxies."%(ndisc,nbulge))

#edisc = imk.e_disc(ndisc,seed=seed_edisc,figout=figdir+keyn+"_eDisc.png")
#ebulge = imk.e_bulge(nbulge, seed=seed_ebulge,figout=figdir+keyn+"_eBulge.png")
edisc = imk.e_disc(ndisc,seed=seed_edisc,figout=None)
ebulge = imk.e_bulge(nbulge, seed=seed_ebulge,figout=None)
#hlrdisc = 1.67835*imk.re_disk(ndisc,mag_gal,rlim=[0.1,2.7],seed=seed_hdisc,figout=figdir+keyn+"_scaleLen.png")
#hlrbulge = imk.re_disk(nbulge,mag_gal,rlim=[0.1,2.7],seed=seed_hbulge)

hlrdisc = imk.re_disk2(mag_gal[id0]) # actually this is scalelength
hlrbulge = imk.re_disk2(mag_gal[id1])/0.6 # equivalent to the scalelength of disc galaxies

logger.info("    The figures of these prior distributions are saved into %s."%params["dir_fig"])
logger.info("")

ell[id0]=edisc; ell[id1]=ebulge
theta[:ngal] = imk.theta(ngal,seed=seed_theta)
hlr[id0]=hlrdisc; hlr[id1]=hlrbulge
bfrac[:ngal] = bfracGal

for i in range(nobj):
    gal_shape = galsim.Shear(g=ell[i], beta=theta[i]*galsim.degrees)
    ell1[i], ell2[i] = gal_shape.g1, gal_shape.g2
    if shear_method == "ps_shear":
        wp = galsim.CelestialCoord(ra[i]*galsim.degrees,dec[i]*galsim.degrees)
        image_posi = wcs.toImage(wp)
        world_posi = affine.toWorld(image_posi)
        g1[i], g2[i], mui = ps.getLensing(pos = world_posi)
        gtt = galsim.Shear(g1=g1[i],g2=g2[i])
        gt[i], beta[i] = gtt.g, gtt.beta.rad()*180.0/np.pi
    elif shear_method == 'constant':
        gtt = galsim.Shear(g1=g1[i],g2=g2[i])
        gt[i], beta[i] = gtt.g, gtt.beta.rad()*180.0/np.pi
        #gtt = galsim.Shear(g=gt[i],beta=beta[i]*galsim.degrees)
        #g1[i], g2[i] = gtt.g1, gtt.g2
    else: # shear_method=="extra":
        gtt = galsim.Shear(g1=g1[i],g2=g2[i])
        gt[i], beta[i] = gtt.g, gtt.beta.rad()*180.0/np.pi

#ell1 += -0.05
#ell2 += 0.01
# observed ellipticity
eg1[:ngal],eg2[:ngal],egt[:ngal],egd[:ngal] = imk.eObs(ell1[:ngal],ell2[:ngal],g1[:ngal],g2[:ngal])

e1_mean = np.mean(ell1[:ngal])
e2_mean = np.mean(ell2[:ngal])
e1_std  = np.std(ell1[:ngal])
e2_std  = np.std(ell2[:ngal])
g1_mean = np.mean(g1[:ngal])
g2_mean = np.mean(g2[:ngal])
g1_std  = np.std(g1[:ngal])
g2_std  = np.std(g2[:ngal])

eg1_mean = np.mean(eg1[:ngal])
eg2_mean = np.mean(eg2[:ngal])
eg1_std  = np.std(eg1[:ngal])
eg2_std  = np.std(eg2[:ngal])

logger.info("    Mean (e1,e2) = (%.6f, %.6f);"%(e1_mean, e2_mean))
logger.info("    SD   (e1,e2) = (%.6f, %.6f);"%(e1_std,  e2_std))
logger.info("    Mean (g1,g2) = (%.6f, %.6f)."%(g1_mean, g2_mean))
logger.info("    SD   (g1,g2) = (%.6f, %.6f);"%(g1_std,  g2_std))
logger.info("    Mean (e1_obs,e2_obs) = (%.6f, %.6f)."%(eg1_mean, eg2_mean))
logger.info("    SD   (e1_obs,e2_obs) = (%.6f, %.6f);"%(eg1_std,  eg2_std))
logger.info("    ** NOTE: change the seed parameters if the result is not satisfactory.")
logger.info("")


# write the catalog
table_array = np.array([index,ra,dec,mag,ell,theta,ell1,ell2,gt,beta,g1,g2,egt,egd,eg1,eg2,hlr,bfrac])
table_array = np.transpose(table_array)

logger.info("^_^ Write mock catalog into %s."%params["mock_catalog"]%field)
gsf.write_FitsTable(contn, table_array, outcat=mockn)

t1 = time.time()
tt = (t1-t0)/60.0
logger.info("    DONE in %.2f minutes!"%tt)
confnNew = catdir+confn+"_%s"%field
os.popen("cp %s %s"%(confnAbs,confnNew))
logger.info("^_^ The configuration file and simulated images are saved into %s."%(confn+"_%s"%field))
logger.info("^_^ All files are saved to directory %s."%params["dir_cat"])
