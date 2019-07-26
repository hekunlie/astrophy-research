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


omega_m0 = float(argv[1])
omega_lambda0 = 1 - omega_m0
h = float(argv[2])
z = float(argv[3])
ang = float(argv[4])# arcmin

radian = ang/60*numpy.pi/180

C_0_hat = 2.99792458
H_0 = 100*h
coeff = 1000*C_0_hat/h

cosmos = FlatLambdaCDM(H_0, Om0=omega_m0)

com_dist = cosmos.comoving_distance(z).value
scale_len = radian*com_dist

print("Z: %.5f. H0: %.3f. Omega_m: %.3f\nCom_dist: %f Mpc/h ( %f Mpc)."%(z, H_0, omega_m0,com_dist*h,com_dist))
print("%.3f arcmin (%.5f rad) ~ %.5f Mpc"%(ang, radian, scale_len))