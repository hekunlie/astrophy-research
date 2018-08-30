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

with open("%s/work/envs/envs.dat"%my_home, "r") as f:
    contents = f.readlines()
for path in contents:
    if "=" in path:
        env_location, env_path = path.split("=")[0:2]
        if "cf_data_path" == env_location:
            cf_data_path = env_path


stamp_size = 64
pts_num = 45
radius = 8
psf_r = 4
flux = 100

mid_rank = int(cpus/2)
num_total = 1000000
num_part = int(num_total/(cpus/2))

if rank < mid_rank:
    g = numpy.load("./cor_g.npz")['arr_0'][:, :2]
    ellip_total = numpy.load("./ellp.npz")["arr_0"]
else:
    g = numpy.load("./cor_g.npz")['arr_0'][:, 2:4]
    ellip_total = numpy.load("./ellp.npz")["arr_1"]

g1 = g[:, 0][rank*num_part:(rank+1)*num_part]
g2 = g[:, 1][rank*num_part:(rank+1)*num_part]

e1 = ellip_total[:, 0][rank*num_part:(rank+1)*num_part]
e2 = ellip_total[:, 1][rank*num_part:(rank+1)*num_part]

data = numpy.zeros((num_part, 4))
fq = Fourier_Quad(stamp_size, rank*10+123)
psf_img = fq.cre_psf(4, 1, "Moffat")
gal_pool = []

for i in range(num_part):
    pts = fq.ran_pos(45, 8, (g1[i], g2[i]))
    gal_img = fq.convolve_psf(pts, 4, 100, "Moffat")
    res = fq.shear_est(gal_img, psf_img)
    data[i] = res[0], res[1], res[2], res[3]

    if rank == 0 and i < 100:
        gal_pool.append(gal_img)

if rank == 0:
    imgs = fq.stack(gal_pool,10)
    hdu = fits.PrimaryHDU(imgs)
    hdu.writeto("./gal.fits",overwrite=True)


if rank == 0:
    recv_buffer = numpy.empty((2*num_total,4))
else:
    recv_buffer = None
comm.Gather(data,recv_buffer, root=0)

if rankk == 0:

    data_name = cf_data_path + "pts.npz"
    numpy.savez(cf_data_path, recv_buffer)

    pic_name = cf_data_path + "pts_fig.png"
    plt.figure(figsize=(8,8))
    plt.hist2d(recv_buffer[:num_total,0], recv_buffer[num_total:num_total*2,0], 5)
    plt.savefig(pic_name)
    plt.close()






