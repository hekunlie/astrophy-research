import matplotlib
matplotlib.use('Agg')
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path
path.append('%s/work/fourier_quad/'%my_home)
import time
from mpi4py import MPI
from Fourier_Quad import *
import tool_box
import numpy
import h5py
import warnings
from sys import argv
import galsim
from astropy.io import fits
import matplotlib.pyplot as plt

cmd = argv[1]

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

t1 = time.time()

envs_path = "%s/work/envs/envs.dat"%my_home
data_path, simu_para_path = tool_box.config(envs_path, ['get','get'], [["correlation", "simu_path", "1"],
                                                                       ["correlation", "simu_para_path", "1"]])
para_file = simu_para_path + "para.ini"
para_items = tool_box.config(para_file, ['get', 'get','get','get'], [["para", 'num', 1], ['para', 'size', 1],
                                                                     ["para", 'cov11', 1], ['para', 'cov12', 1]])
total_num, stamp_size = int(para_items[0]), int(para_items[1])
cov11, cov12 = float(para_items[2]), float(para_items[3])

cpu_mid = int(cpus/2)
region, data_tag = divmod(rank, cpu_mid)
partial_num = int(total_num/cpu_mid)
start_pt, end_pt = partial_num*rank, partial_num*(rank + 1)

para_h5path = simu_para_path + "para.hdf5"

if cmd == "para":
    if rank == 0:
        para_h5 = h5py.File(para_h5path, "w")

        for i in range(2):
            cov = [[cov11, cov12+i*0.0002], [cov12+i*0.0002, cov11]]
            gs = tool_box.rand_gauss2n(num=total_num, means=[0, 0], cov=cov).T
            print(gs.shape)
            e1, e2 = tool_box.ellip_mock(total_num, 10*i+1)
            plt.subplot(121)
            plt.hist(e1, 100)
            plt.subplot(122)
            plt.hist(e2, 100)
            plt.savefig(data_path+"ellip_%d.png"%i)
            plt.close()
            plt.scatter(gs[:,0], gs[:,1], s=1)
            plt.xlim(-0.1, 0.1)
            plt.ylim(-0.1, 0.1)
            plt.savefig(data_path+"gg_%d.png"%i)
            plt.close()
            for j in range(cpu_mid):
                g_set_1 = "/0/g%d/%d/"%(i+1, j)
                g_set_2 = "/1/g%d/%d/"%(i+1, j)
                para_h5[g_set_1] = gs[partial_num*j: partial_num*(j+1), 0]
                para_h5[g_set_2] = gs[partial_num*j: partial_num*(j+1), 1]
                e_set_1 = "/%d/e1/%d/"%(i, j)
                e_set_2 = "/%d/e2/%d/"%(i, j)
                para_h5[e_set_1] = e1[partial_num*j: partial_num*(j+1)]
                para_h5[e_set_2] = e2[partial_num*j: partial_num*(j+1)]
        para_h5.close()

else:
    para_h5 = h5py.File(para_h5path,"r")

    g1_setname = "/%d/g1/%d/"%(region, data_tag)
    g1 = para_h5[g1_setname].value
    g2_setname = "/%d/g2/%d/"%(region, data_tag)
    g2 = para_h5[g2_setname].value
    e1_setname = "/%d/e1/%d/"%(region, data_tag)
    e1 = para_h5[e1_setname].value
    e2_setname = "/%d/e2/%d/"%(region, data_tag)
    e2 = para_h5[e2_setname].value
    para_h5.close()

    fq = Fourier_Quad(stamp_size, rank*10+10)

    gal_radius = 1.0
    pixel_scale = 0.2
    gal_flux = 1.e5

    psf = galsim.Moffat(beta=3.5, scale_radius=0.8, flux=1.0, trunc=3)
    psf_img = galsim.ImageD(stamp_size, stamp_size)
    psf.drawImage(image=psf_img, scale=pixel_scale)
    psf_pow = fq.pow_spec(psf_img.array)
    hlr = fq.get_radius_new(psf_pow, 2)

    gal_pool = []
    sp = (partial_num, 4)
    measure_data = numpy.zeros(sp)
    for i in range(partial_num):

        gal = galsim.Sersic(half_light_radius=gal_radius, n=2, trunc=4.5 * gal_radius)  # be careful
        gal = gal.shear(e1=e1[i], e2=e2[i])

        gal_s = gal.withFlux(gal_flux)
        gal_g = gal_s.shear(g1=g1[i], g2=g2[i])
        gal_c = galsim.Convolve([gal_g, psf])
        img = galsim.ImageD(stamp_size, stamp_size)
        gal_c.drawImage(image=img, scale=pixel_scale)
        gal_img = img.array# + fq.draw_noise(0, noise_sig)
        if rank == 0 and i < 10000:
            gal_pool.append(gal_img)

        measure_data[i] = fq.shear_est(gal_img, psf_pow, F=True)[:4]

    rev_buffer = None
    if rank == 0:
        rev_buffer = numpy.empty((2*total_num, 4), dtype="numpy.float64")
    comm.Gather(measure_data, rev_buffer, root=0)

    if rank == 0:
        galaxies = fq.stack(gal_pool, 100)
        hdu = fits.PrimaryHDU(galaxies)
        hdu.writeto(data_path + "gal.fits", overwrite=True)
        data_1, data_2 = rev_buffer[:total_num], rev_buffer[total_num:]
        numpy.savez(data_path+"data.npz", data_1, data_2)
t2 = time.time()
if rank == 0:
    print(t2 - t1)