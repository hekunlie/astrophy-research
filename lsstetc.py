"""An exposure time calculator for LSST.  Uses GalSim to draw a galaxy with specified magnitude,
shape, etc, and then uses the same image as the optimal weight function.  Derived from D. Kirkby's
notes on deblending.
"""

import numpy as np

#import galsim

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

    # def draw(self, profile, add_noise=None):
    #     img = galsim.ImageD(self.stamp_size, self.stamp_size, scale=self.pixel_scale)
    #     #profile = profile.withFlux(self.flux(mag))
    #     profile.drawImage(image=img)
    #     image = img.array
    #     signal = np.sum(image**2)
    #     noise = np.sqrt(signal * self.sky)
    #     snr = signal/noise
    #     if add_noise is not None:
    #         dirt_img = image + np.random.normal(loc=0, scale=self.sigma_sky, size=self.stamp_size**2).reshape(self.stamp_size, self.stamp_size)
    #         return dirt_img, snr
    #     else:
    #         return image

    # def SNR(self, profile, mag):
    #     img = self.draw(profile, add_noise=None)
    #     mask = img.array > (self.threshold * self.sigma_sky)
    #     imgsqr = img.array**2*mask
    #     signal = imgsqr.sum()
    #     noise = np.sqrt((imgsqr * self.sky).sum())
    #     return signal / noise
    #
    # def err(self, profile, mag):
    #     snr = self.SNR(profile, mag)
    #     return 2.5 / np.log(10) / snr
    #
    # def display(self, profile, mag, noise=True):
    #     img = self.draw(profile, noise)
    #     import matplotlib.pyplot as plt
    #     import matplotlib.cm as cm
    #     plt.imshow(img.array, cmap=cm.Greens)
    #     plt.colorbar()
    #     plt.show()




