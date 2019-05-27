import astropy.wcs as wcs
import astropy.units as U

a = fits.open("...")
awcs = wcs.WCS(a[0].header)
cx, cy = awcs.all_world2pix(ra*U.deg, dec*U.deg, 0)
