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
import galsim
import matplotlib.pyplot as plt
from astropy.io import fits
import tool_box
from mpi4py import MPI
import h5py
import configparser


comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

config = configparser.ConfigParser()
config.read("%s/work/envs/envs.dat"%my_home)

cf_data_path = config.get("cfht", "data_path_simu")

time_s = time.time()

with open("./paras.dat","r") as f:
    contents = f.readlines()
for item in contents:
    if ":" in item:
        item_name, item_val = item.split(":")
        if "stamp_size" == item_name:
            stamp_size = int(item_val)
        if "pts_num" == item_name:
            pts_num = int(item_val)
        if "max_radius" == item_name:
            radius = float(item_val)
        if "psf_r" == item_name:
            psf_r = float(item_val)
        if "total_num" == item_name:
            num_total = int(item_val)

flux = 100

mid_rank = int(cpus/2.)
num_part = int(num_total/(cpus/2.))

para_path = cf_data_path + "para.hdf5"
time.sleep(rank*1.5)
paras = h5py.File(para_path, "r")
ellips = paras["/ellips"].value
# mags = paras["/magnitudes"].value
shear = paras["/shear"].value

if rank < mid_rank:
    total_g = shear[:,0]
    ellip_total = ellips[: num_total]
    # mag_total = mags[:num_total]

else:
    total_g = shear[:,1]
    ellip_total = ellips[num_total: 2*num_total]
    # mag_total = mags[num_total: 2*num_total]

rank_tag = divmod(rank, mid_rank)[1]
tag_s, tag_e = rank_tag*num_part, (rank_tag+1)*num_part

g1 = total_g[0:num_total][tag_s:tag_e]
g2 = total_g[num_total: 2*num_total][tag_s:tag_e]

e1 = ellip_total[tag_s:tag_e, 0]
e2 = ellip_total[tag_s:tag_e, 1]
# mag = mag_total[tag_s:tag_e]

data = numpy.zeros((num_part, 4))
fq = Fourier_Quad(stamp_size, rank*10+1234)
# psf_img = fq.cre_psf(4, 1, "Moffat")
gal_pool = []

print(rank, tag_s, tag_e,g1.shape, g2.shape)
pixel_scale = 0.2
psf = galsim.Moffat(beta=3.5, scale_radius=1.0, flux=1.0, trunc=4.5)
psf_ = galsim.ImageD(stamp_size, stamp_size)
psf.drawImage(image=psf_, scale=pixel_scale)
psf_img = psf_.array
# if rank < mid_rank: for testing
#     ig1, ig2 = 0.02, -0.03
# else:
#     ig1, ig2 = -0.01, 0.023

for i in range(num_part):

    # pts = fq.ran_pos(45, 8, (g1[i], g2[i]))[1]
    # gal_img = fq.convolve_psf(pts, 4, 100, "Moffat")

    ra = numpy.random.uniform(0.8, 1.4, 1)[0]
    alp = numpy.random.uniform(0.01,2,1)[0]
    bulge = galsim.Sersic(half_light_radius=ra, n=3.5, trunc=4.5 * ra)  # be careful
    # disk = galsim.Sersic(half_light_radius=ra, n=1, trunc=4.5 * ra)  # be careful
    # gal = bulge * btr + disk * (1 - btr)
    # gal = gal.shear(e1=e1, e2=e2)
    gal = bulge.shear(e1=e1[i], e2=e2[i])
    gal_s = gal.withFlux(1.e6*alp)
    gal_g = gal_s.shear(g1=g1[i], g2=g2[i])
    gal_c = galsim.Convolve([gal_g, psf])
    img = galsim.ImageD(stamp_size, stamp_size)
    gal_c.drawImage(image=img, scale=pixel_scale)
    gal_img = img.array# + fq.draw_noise(0, noise_sig)

    res = fq.shear_est(gal_img, psf_img)
    data[i] = res[0], res[1], res[2], res[3]

    if rank == 0 and i < 100:
        gal_pool.append(gal_img)

if rank == 0:
    imgs = fq.stack(gal_pool,10)
    hdu = fits.PrimaryHDU(imgs)
    hdu.writeto("./gal.fits",overwrite=True)


if rank == 0:
    recv_buffer = numpy.empty((2*num_total, 4))
else:
    recv_buffer = None
comm.Gather(data,recv_buffer, root=0)

time_e = time.time()
if rank == 0:

    print("Time used: ", time_e - time_s)
    data_name = cf_data_path + "gal.npz"
    numpy.savez(data_name, recv_buffer)

    pic_name = cf_data_path + "gal_fig.png"
    plt.figure(figsize=(8,8))
    plt.hist2d(recv_buffer[:num_total,0], recv_buffer[num_total:num_total*2,0], 100)
    plt.savefig(pic_name)
    plt.close()






