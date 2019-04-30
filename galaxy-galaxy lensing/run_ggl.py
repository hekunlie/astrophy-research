from subprocess import Popen
import time
from sys import argv
import os
from mpi4py import MPI



comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

area = [1,3,4]
radius_label = [i for i in range(13)]

for area_id in area:
    cmd = "mpirun -n 1 ./ggl_test %d %d"%(area_id, rank)
    a = Popen(cmd, shell=True)
    a.wait()