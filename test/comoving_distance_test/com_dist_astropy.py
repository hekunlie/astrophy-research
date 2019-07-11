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
h = float(argv[3])

C_0_hat = 2.99792458
H_0 = 100*h
coeff = 1000*C_0_hat/h

num = 1000

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)

redshift = numpy.linspace(0, 10, 1000).reshape(num, 1)
com_dist = numpy.zeros((num, 1))

if cmd == "generate":
    for i in range(num):
        com_dist[i,0] = cosmos.comoving_distance(redshift[i]).value

    h5f = h5py.File("dist.hdf5","w")
    h5f["/Z"] = redshift
    h5f["/DISTANCE"] = com_dist
    h5f["/PARAMETER"] = numpy.array([omega_m0, omega_lambda0, h, C_0_hat])
    h5f.close()

    cmd = "./dist_test"
    a = Popen(cmd, shell=True)
    a.wait()

if cmd == "compare":

    h5f = h5py.File("dist.hdf5","r")
    dist_astropy = h5f["/DISTANCE"].value[:, 0]
    dist_cpp = h5f["/DISTANCE_C"].value[:, 0]
    h5f.close()

    img = Image_Plot()
    img.subplots(1, 2)

    img.axs[0][0].scatter(redshift, dist_astropy, label="TRUE by Astropy")
    img.axs[0][0].scatter(redshift, dist_cpp, label="CPP")
    img.axs[0][1].hist(dist_astropy-dist_cpp, 100)
    img.axs[0][0].legend()
    img.set_label(0,0,0,"Com_Dist")
    img.set_label(0,0,1,"Z")
    img.set_label(0,1,0,"Number")
    img.set_label(0,1,1,"Difference")
    img.save_img("compare_result.png")


