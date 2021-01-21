import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
from subprocess import Popen
import numpy
from astropy.coordinates import SkyCoord
from astropy import units as astro_unit


for i in range(10):
    ra1, ra2 = numpy.random.uniform(0,360,2).tolist()
    dec1, dec2 = numpy.random.uniform(-90,90,2).tolist()

    c1 = SkyCoord(ra=ra1 * astro_unit.deg, dec=dec1 * astro_unit.deg, frame='fk5')
    c2 = SkyCoord(ra=ra2 * astro_unit.deg, dec=dec2 * astro_unit.deg, frame='fk5')
    sep_angle = c1.separation(c2).radian
    print(sep_angle)

    cmd = "./separation_angle_test %f %f %f %f"%(ra1, dec1, ra2, dec2)
    a = Popen(cmd, shell=True)
    a.wait()
    print("\n")
