from astropy.coordinates import SkyCoord
from astropy import units
from sys import argv

ra1, dec1 = float(argv[1]), float(argv[2])
ra2, dec2 = float(argv[3]), float(argv[4])

c1 = SkyCoord(ra=ra1*units.deg, dec=dec1*units.deg,frame="fk5")
c2 = SkyCoord(ra=ra2*units.deg, dec=dec2*units.deg,frame="fk5")
sep = c1.separation(c2)

print("(%10.5f,%10.5f) <-- %10.5f rad (%10.5f deg) --> (%10.5f,%10.5f)"%(ra1, dec1,sep.radian, sep.deg,ra2, dec2))
