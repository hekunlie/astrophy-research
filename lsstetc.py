"""An exposure time calculator for LSST.  Uses GalSim to draw a galaxy with specified magnitude,
shape, etc, and then uses the same image as the optimal weight function.  Derived from D. Kirkby's
notes on deblending.
"""

import numpy as np

import galsim

# Some constants
# --------------
#
# LSST effective area in meters^2
A = 319/9.6  # etendue / FoV.  I *think* this includes vignetting

# zeropoints from DK notes in photons per second per pixel
# should eventually compute these on the fly from filter throughput functions.
s0 = {'u': A*0.732,
      'g': A*2.124,
      'r': A*1.681,
      'i': A*1.249,
      'z': A*0.862,
      'Y': A*0.452}
# Sky brightnesses in AB mag / arcsec^2.
# stole these from http://www.lsst.org/files/docs/gee_137.28.pdf
# should eventually construct a sky SED (varies with the moon phase) and integrate to get these
B = {'u': 22.8,
     'g': 22.2,
     'r': 21.3,
     'i': 20.3,
     'z': 19.1,
     'Y': 18.1}
# number of visits
# From LSST Science Book
fiducial_nvisits = {'u': 56,
                    'g': 80,
                    'r': 180,
                    'i': 180,
                    'z': 164,
                    'Y': 164}
# exposure time per visit
visit_time = 30.0
# Sky brightness per arcsec^2 per second
sbar = {}
for k in B:
    sbar[k] = s0[k] * 10**(-0.4*(B[k]-24.0))

# And some random numbers for drawing



class ETC(object):
    def __init__(self, band, pixel_scale=None, stamp_size=None, threshold=0.0,
                 nvisits=None):
        self.pixel_scale = pixel_scale
        self.stamp_size = stamp_size
        self.threshold = threshold
        self.band = band
        if nvisits is None:
            self.exptime = fiducial_nvisits[band] * visit_time
        else:
            self.exptime = nvisits * visit_time
        self.sky = sbar[band] * self.exptime * self.pixel_scale**2
        self.sigma_sky = np.sqrt(self.sky)
        self.s0 = s0[band]

    def flux(self, mag):
        return self.s0 * 10**(-0.4*(mag - 24.0)) * self.exptime

    def draw(self, profile, mag, seed=None, add_noise=None):
        img = galsim.ImageD(self.stamp_size, self.stamp_size, scale=self.pixel_scale)
        profile = profile.withFlux(self.flux(mag))
        profile.drawImage(image=img)
        signal = np.sum(img.array**2)
        noise = np.sqrt((signal * self.sky).sum())
        snr = signal/noise
        if add_noise is not None:
            bd = galsim.BaseDeviate(seed)
            gd = galsim.GaussianNoise(bd, sigma = self.sigma_sky)
            img.addNoise(gd)
            noise_img = galsim.ImageD(self.stamp_size, self.stamp_size, scale=self.pixel_scale)
            bdn = galsim.BaseDeviate(seed+2)
            gdn = galsim.GaussianNoise(bdn, sigma = self.sigma_sky)
            noise_img.addNoise(gdn)
            return img, noise_img,snr
        else:
            return img

    def SNR(self, profile, mag):
        img = self.draw(profile, mag,noise=False)
        mask = img.array > (self.threshold * self.sigma_sky)
        imgsqr = img.array**2*mask
        signal = imgsqr.sum()
        noise = np.sqrt((imgsqr * self.sky).sum())
        return signal / noise

    def err(self, profile, mag):
        snr = self.SNR(profile, mag)
        return 2.5 / np.log(10) / snr

    def display(self, profile, mag, noise=True):
        img = self.draw(profile, mag, noise)
        import matplotlib.pyplot as plt
        import matplotlib.cm as cm
        plt.imshow(img.array, cmap=cm.Greens)
        plt.colorbar()
        plt.show()


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser()

    # Filter
    parser.add_argument("--band", default='i',
                        help="band for simulation (Default 'i')")

    # PSF structural arguments
    PSF_profile = parser.add_mutually_exclusive_group()
    PSF_profile.add_argument("--kolmogorov", action="store_true",
                             help="Use Kolmogorov PSF (Default Gaussian)")
    PSF_profile.add_argument("--moffat", action="store_true",
                             help="Use Moffat PSF (Default Gaussian)")
    parser.add_argument("--PSF_beta", type=float, default=3.0,
                        help="Set beta parameter of Moffat profile PSF. (Default 2.5)")
    parser.add_argument("--PSF_FWHM", type=float, default=0.67,
                        help="Set FWHM of PSF in arcsec (Default 0.67).")
    parser.add_argument("--PSF_phi", type=float, default=0.0,
                        help="Set position angle of PSF in degrees (Default 0.0).")
    parser.add_argument("--PSF_ellip", type=float, default=0.0,
                        help="Set ellipticity of PSF (Default 0.0)")

    # Galaxy structural arguments
    parser.add_argument("-n", "--sersic_n", type=float, default=1.0,
                        help="Sersic index (Default 1.0)")
    parser.add_argument("--gal_ellip", type=float, default=0.3,
                        help="Set ellipticity of galaxy (Default 0.3)")
    parser.add_argument("--gal_phi", type=float, default=0.0,
                        help="Set position angle of galaxy in radians (Default 0.0)")
    parser.add_argument("--gal_HLR", type=float, default=0.2,
                        help="Set galaxy half-light-radius. (default 0.5 arcsec)")

    # Simulation input arguments
    parser.add_argument("--pixel_scale", type=float, default=0.2,
                        help="Set pixel scale in arcseconds (Default 0.2)")
    parser.add_argument("--stamp_size", type=int, default=31,
                        help="Set postage stamp size in pixels (Default 31)")

    # Magnitude!
    parser.add_argument("--mag", type=float, default=25.3,
                        help="magnitude of galaxy")
    # threshold
    parser.add_argument("--threshold", type=float, default=0.0,
                        help="Threshold, in sigma-sky units, above which to include pixels")

    # Observation characteristics
    parser.add_argument("--nvisits", type=int, default=None)

    # draw the image!
    parser.add_argument("--display", action='store_true',
                        help="Display image used to compute SNR.")

    args = parser.parse_args()

    if args.kolmogorov:
        psf = galsim.Kolmogorov(fwhm=args.PSF_FWHM)
    elif args.moffat:
        psf = galsim.Moffat(fwhm=args.PSF_FWHM, beta=args.PSF_beta)
    else:
        psf = galsim.Gaussian(fwhm=args.PSF_FWHM)
    psf = psf.shear(e=args.PSF_ellip, beta=args.PSF_phi*galsim.radians)

    gal = galsim.Sersic(n=args.sersic_n, half_light_radius=args.gal_HLR)
    gal = gal.shear(e=args.gal_ellip, beta=args.gal_phi*galsim.radians)

    profile = galsim.Convolve(psf, gal)

    etc = ETC(args.band, pixel_scale=args.pixel_scale, stamp_size=args.stamp_size,
              threshold=args.threshold, nvisits=args.nvisits)

    print
    print "input"
    print "------"
    print "band: {}".format(args.band)
    print "magnitude: {}".format(args.mag)
    print
    print "output"
    print "------"
    print "signal: {}".format(etc.flux(args.mag))
    print "sky noise: {}".format(etc.sigma_sky)
    print "SNR: {}".format(etc.SNR(profile, args.mag))
    print "mag err: {}".format(etc.err(profile, args.mag))

    if args.display:
        etc.display(profile, args.mag)
