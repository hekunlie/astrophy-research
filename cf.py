import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
# path.append("E:/Github/astrophy-research/")
import time
from Fourier_Quad import Fourier_Quad
# import galsim
import matplotlib.pyplot as plt
from astropy.io import fits
import tool_box
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

if rank < 10:
    g = numpy.load("./cor_g.npz")['arr_0'][:, :2]
else:
    g = numpy.load("./cor_g.npz")['arr_0'][:, 2:4]

g1 = g[:, 0]
g2 = g[:, 1]

