import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
from mpi4py import MPI
import h5py
from Fourier_Quad import Fourier_Quad
import time


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()

t1 = time.time()

stamp_size = 44
pts_num = 30
max_radius = 7
psf_scale = 4
psf_flux = 1
gal_flux = 5000
noise_sig = 60
g1 = -0.04
g2 = 0.04
psf_type = "Moffat"

total_num = 400000

step_size = int(argv[1])
seed = int(argv[2])


itemsize = MPI.DOUBLE.Get_size()

if rank == 0:

    nbytes = numprocs*total_num*10*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
data = numpy.ndarray(buffer=buf1, dtype='d', shape=(total_num*numprocs, 10)) # array filled with zero


fq = Fourier_Quad(stamp_size, seed)

psf_img = fq.cre_psf(psf_scale, psf_flux, psf_type)
psf_pow = fq.pow_spec(psf_img)
fq.get_radius_new(psf_pow, 2)

pts = fq.ran_pts(pts_num, max_radius, step_size)
pst_s = fq.shear(pts, g1, g2)

gal_img_nf = fq.convolve_psf(pst_s, psf_scale, gal_flux, psf_type)
gal_pow_nf = fq.pow_spec(gal_img_nf)

for i in range(total_num):

    noise_1 = fq.draw_noise(0, noise_sig)
    noise_pow = fq.pow_spec(noise_1)

    gal_img_n = gal_img_nf + noise_1
    gal_pow_n = fq.pow_spec(gal_img_n)

    cross_term = gal_pow_n - gal_pow_nf - noise_pow

    data[rank*total_num + i,:5] = fq.shear_est(gal_pow_nf, psf_pow, F=True)
    data[rank*total_num + i,5:10] = fq.shear_est(cross_term, psf_pow, F=True)

comm.Barrier()
t2 = time.time()

if rank == 0:

    data_nf = numpy.zeros((1,5))
    data_nf[0] = fq.shear_est(gal_pow_nf, psf_pow, F=True)

    h5f = h5py.File("./CT_data_%.2f.hdf5"%step_size, "w")
    h5f["/data"] = data
    h5f["/data_noise_free"] = data_nf
    h5f.close()

    print(t2 - t1)