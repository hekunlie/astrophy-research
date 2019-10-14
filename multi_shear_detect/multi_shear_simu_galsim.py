import numpy
from sys import path,argv
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import galsim
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

total_num = int(argv[1])

g1 = [-0.04, 0.01, 0.03]
g2 = [0.04, -0.01, -0.03]

shear_num = len(g1)

