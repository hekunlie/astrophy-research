import numpy
from astropy.coordinates import SkyCoord
from astropy import units
import h5py

num = 20000
num2 = int(num/2)
ra = numpy.random.uniform(0, 360, num).reshape((num2,2))
dec = numpy.random.uniform(-90, 90, num).reshape((num2,2))

pts = numpy.zeros((num2, 6))
pts[:,:2] = ra[:,:2]
pts[:,2:4] = dec[:,:2]

for i in range(num2):
    c1 = SkyCoord(ra=pts[i,0]*units.deg,dec=pts[i,2]*units.deg,frame ="fk5")
    c2 = SkyCoord(ra=pts[i,1]*units.deg,dec=pts[i,3]*units.deg,frame ="fk5")
    sep = c1.separation(c2)
    pts[i,4] = sep.radian
    pts[i,5] = sep.deg
    print("(%10.5f,%10.5f) <-- %10.5f rad (%10.5f deg) --> (%10.5f,%10.5f)"%(pts[i,0],pts[i,2],pts[i,4],pts[i,5],pts[i,1],pts[i,3]))

h5f = h5py.File("Sep_test.hdf5","w")
h5f["/data"] = pts
h5f.close()


