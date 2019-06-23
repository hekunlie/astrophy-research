import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/'%my_home)
import tool_box
from astropy.cosmology import FlatLambdaCDM
import h5py


# input the Z value for checking
z_check = float(argv[1])

cosmos = FlatLambdaCDM(H0=70, Om0=0.31)
# comoving distance
com_dist = cosmos.comoving_distance(z_check).value

H0 = 70
C = 299792.458
coeff = C/H0
# the comoving distance (only the integrate part) calculated with Omega_M = 0.31, Omega_Lanbda = 0.69
f = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5","r")
z = f["/redshift"].value
dist = f["/distance"].value
f.close()

tag = tool_box.find_near(z, z_check)
co_distance = dist[tag]*coeff
print("Z: %7.5f, Dist: %12.6f Mpc, Dist (astropy): %12.6f Mpc. Diff: %f Mpc"%(z[tag], co_distance, com_dist, co_distance - com_dist))