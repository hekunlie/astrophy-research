import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import numpy
from astropy.coordinates import SkyCoord
from astropy import units
import h5py
from sys import path, argv
path.append("E:/Github/astrophy-research/mylib/")
path.append('%s/work/mylib/' % my_home)
from plot_tool import Image_Plot
from subprocess import Popen


# generate 20000 galaxy pairs and calculate the separation angle
# for separation() debugging
cmd = argv[1]

if cmd == "generate":
    num = 10000
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

    h5f = h5py.File("./data/sep_test.hdf5","w")
    h5f["/data"] = pts
    h5f.close()

    cmd = "./FQ_sep_test 1 1 0 0 1"
    a = Popen(cmd, shell=True)
    a.wait()


if cmd == "compare":

    f = h5py.File("./data/sep_test_cpp.hdf5", "r")
    data1 = f["/data"].value
    f.close()

    f = h5py.File("./data/sep_test.hdf5", "r")
    data2 = f["/data"].value
    f.close()

    img = Image_Plot()
    img.subplots(1, 2)

    img.axs[0][0].hist(data1[:,0]-data2[:,4], 100)
    img.axs[0][1].hist(data1[:,1]-data2[:,5], 100)

    img.save_img("compare_result.png")
    # img.show_img()

