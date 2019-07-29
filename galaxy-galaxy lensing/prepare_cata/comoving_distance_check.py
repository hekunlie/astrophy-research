import matplotlib
matplotlib.use("Agg")
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path, argv
path.append("E:/Github/astrophy-research/mylib/")
path.append('%s/work/mylib/' % my_home)
from plot_tool import Image_Plot
import tool_box
from subprocess import Popen
import numpy
from astropy.cosmology import FlatLambdaCDM
import h5py
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


Omg_m0 = float(argv[1])
Omg_lam0 = 1 - Omg_m0
data_path = argv[2]
set_name = argv[3]
ch_num = int(argv[4])

h5f = h5py.File(data_path,"r")
redshift = h5f["%sZ"%set_name].value

total_num = redshift.shape[0]
# my_list = numpy.array(tool_box.allot([i for i in range(total_num)], cpus)[rank])
my_list = [i for i in range(total_num)]
ch_pts = numpy.random.choice(my_list, ch_num, replace=False).tolist()

com_dist_file = h5f["%sDISTANCE"%set_name].value[ch_pts]
com_dist_integ_file = h5f["%sDISTANCE_INTEG"%set_name].value[ch_pts]
redshift = redshift[ch_pts]

h5f.close()

h = 0.7
C_0_hat = 2.99792458
H_0 = 100*h
coeff = 1000*C_0_hat

cosmos = FlatLambdaCDM(H_0, Om0=Omg_m0)
delta_dist = numpy.zeros((2, ch_num))
astro_dist = numpy.zeros((2, ch_num))
for i in range(ch_num):
    com_dist = cosmos.comoving_distance(redshift[i]).value*h
    integ = com_dist/coeff
    print("%.5f, %.6f"%(redshift[i], com_dist_file[i] - com_dist))
    # delta_dist[0, i] = numpy.abs(com_dist_file[i] - com_dist)
    # delta_dist[1, i] = numpy.abs(com_dist_integ_file[i] - integ)
    #
    # astro_dist[0, i] = com_dist
    # astro_dist[1, i] = integ
#
# dy = delta_dist[0].max() - delta_dist[0].min()
# plt_dy = (delta_dist[0].min() - dy*0.2, delta_dist[0].max() + dy*0.2)
# img = Image_Plot()
# img.subplots(2,2)
# img.axs[0][0].scatter(redshift, com_dist_file, c="C1",s=2)
# img.axs[0][0].scatter(redshift, astro_dist[0], c="C2",s=2)
# img.axs[1][0].scatter(redshift, delta_dist[0], c="C1",s=2)
# img.axs[1][0].set_ylim(plt_dy)
#
# dy = delta_dist[1].max() - delta_dist[1].min()
# plt_dy = (delta_dist[1].min() - dy*0.2, delta_dist[1].max() + dy*0.2)
# img.axs[0][1].scatter(redshift, com_dist_integ_file, c="C1",s=2)
# img.axs[0][1].scatter(redshift, astro_dist[1], c="C2",s=2)
# img.axs[1][1].scatter(redshift, delta_dist[1], c="C1",s=2)
# img.axs[1][1].set_ylim(plt_dy)
#
# img.save_img("./pic/%d.png"%rank)

# for i in range(cpus):
#     if i == rank:
#         print(rank,my_list.min(), my_list.max(), delta_dist.min(), delta_dist.max())
#     comm.Barrier()


