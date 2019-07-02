import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append("E:/Github/astrophy-research/mylib/")
path.append('%s/work/mylib/' % my_home)
from plot_tool import Image_Plot
from subprocess import Popen
import numpy
from astropy.cosmology import FlatLambdaCDM
import h5py


cmd = argv[1]
omega_m0 = float(argv[2])
omega_lambda0 = 1 - omega_m0

C_0_hat = 2.99792458
h = 0.7
H_0 = 100*h
coeff = 1000*C_0_hat/h

num = 1000

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)

redshift = numpy.linspace(0, 10, 1000).reshape(num, 1)
com_dist = numpy.zeros((num, 1))

if cmd == "generate":
    for i in range(num):
        com_dist[i,0] = cosmos.comoving_distance(redshift[i])

    h5f = h5py.File("dist.h5py","w")
    h5f["/Z"] = redshift
    h5f["/DISTANCE"] = com_dist
    h5f.close()

    cmd = "./dist_test"
    a = Popen(cmd, shell=True)
    a.wait()

if cmd == "compare":

    h5f = h5py.File("dist.h5py","r")
    dist_astropy = h5f["/DISTANCE"].value
    dist_cpp = h5f["/DISTANCE_C"].value
    h5f.close()

    img = Image_Plot()
    img.subplots(1, 2)

    img.axs[0][0].scatter(redshift, dist_astropy[:, 0], label="TRUE")
    img.axs[0][0].scatter(redshift, dist_cpp[:, 0], label="CPP")
    img.axs[0][1].hist(dist_astropy[:,0]-dist_cpp[:,0], 100)
    img.axs[0][0].legend()

    img.save_img("compare_result.png")


