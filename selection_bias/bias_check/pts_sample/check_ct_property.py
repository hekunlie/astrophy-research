import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
from mpi4py import MPI
import h5py
from Fourier_Quad import Fourier_Quad
import time
from plot_tool import Image_Plot

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()



cmd = argv[1]

if cmd == "simu":


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

    total_num = 200000

    step_size = int(argv[2])
    seed = int(argv[3])


    itemsize = MPI.FLOAT.Get_size()

    if rank == 0:
        print(itemsize)
        nbytes = numprocs*total_num*10*itemsize
    else:
        nbytes = 0

    win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
    buf1, itemsize = win1.Shared_query(0)
    data = numpy.ndarray(buffer=buf1, dtype='f', shape=(total_num*numprocs, 10)) # array filled with zero


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

        gal_pow_dn = gal_pow_n - noise_pow

        cross_term = gal_pow_n - gal_pow_nf - noise_pow

        data[rank*total_num + i,:5] = fq.shear_est(gal_pow_dn, psf_pow, F=True)
        data[rank*total_num + i,5:10] = fq.shear_est(cross_term, psf_pow, F=True)

    comm.Barrier()
    t2 = time.time()

    if rank == 0:

        data_nf = numpy.zeros((1,5), dtype=numpy.float32)
        data_nf[0] = fq.shear_est(gal_pow_nf, psf_pow, F=True)

        h5f = h5py.File("./CT_data_%.02f.hdf5"%step_size, "w")
        h5f["/data"] = data
        h5f["/data_noise_free"] = data_nf[0]
        h5f["/img_noisy"] = gal_img_nf + fq.draw_noise(0, noise_sig)
        h5f.close()

        print(t2 - t1)

else:
    fq = Fourier_Quad(12, 1243)

    bin_num = 4
    h5f = h5py.File("./CT_data_%.02f.hdf5"%(rank+1), "r")
    data = h5f["/data"][()]
    h5f.close()

    mg1 = data[:,5]
    mg2 = data[:,6]
    mnu1 = data[:,7] + data[:,8]
    mnu2 = data[:,7] - data[:,8]

    gh = numpy.linspace(-0.1,0.1, 41)
    chisq1 = fq.get_chisq_range(mg1, mnu1, bin_num, gh)[1]
    chisq2 = fq.get_chisq_range(mg2, mnu2, bin_num, gh)[1]

    h5f = h5py.File("./CT_data_%.02f_chisq.hdf5" % (rank + 1), "w")
    h5f["/chisq1"] = chisq1
    h5f["/chisq2"] = chisq2
    h5f["/gh"] = gh
    h5f.close()

    comm.Barrier()

    if rank == 0:
        img = Image_Plot(xpad=0.2, ypad=0.2)
        img.subplots(1, 5)
        for i in range(numprocs):
            h5f = h5py.File("./CT_data_%.02f_chisq_bin%d.hdf5" % (i + 1, bin_num), "r")
            chisq1 = h5f["/chisq1"][()]
            chisq2 = h5f["/chisq2"][()]
            gh = h5f["/gh"][()]
            h5f.close()
            img.axs[0][i].plot(gh, chisq1, label="g1")
            img.axs[0][i].plot(gh, chisq2, label="g2")
            img.axs[0][i].legend()
        img.save_img("./chisq_bin%d.png"%bin_num)

    comm.Barrier()