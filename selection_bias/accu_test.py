from sys import path, argv
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import tool_box
import numpy
import h5py
import galsim
from Fourier_Quad import Fourier_Quad
from mpi4py import MPI
import time

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


psf_tag = argv[1]

shear_num = 15
g = numpy.linspace(-0.04,0.04,shear_num)
theta = numpy.random.uniform(0,numpy.pi*2,shear_num)
g1 = g*numpy.cos(2*theta)
g2 = g*numpy.sin(2*theta)

stamp_size = 96

fq = Fourier_Quad(stamp_size,432)

# parameters
si_num = 11
sr_num = 11
sersic_index = numpy.linspace(0.4, 2, si_num)
scale_radius = numpy.linspace(0.2, 1, sr_num)

pixel_scale = 0.187

psf_e = 0.1
psf_scale = 0.8

gal_num = int(argv[2])
gal_e = numpy.random.uniform(0.1, 0.7, gal_num)
gal_theta = numpy.random.uniform(0, numpy.pi * 2, gal_num)
gal_e1 = gal_e * numpy.cos(2 * gal_theta)
gal_e2 = gal_e * numpy.sin(2 * gal_theta)


itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = 4*si_num*sr_num*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
mcs = numpy.ndarray(buffer=buf1, dtype='d', shape=(4, si_num, sr_num))


# PSF
if psf_tag == "c":
    psf = galsim.Moffat(beta=3.5, fwhm=psf_scale, flux=1.0, trunc=1.4)
    if rank == 0:
        print("Circle PSF")
else:
    psf = galsim.Moffat(beta=3.5, fwhm=psf_scale, flux=1.0, trunc=1.4).shear(e1=psf_e, e2=0)
    if rank == 0:
        print("Ellipticle PSF")
psf_img = galsim.ImageD(stamp_size, stamp_size)
psf.drawImage(image=psf_img, scale=pixel_scale)
psf_img = psf_img.array
psf_pow = fq.pow_spec(psf_img)

fq.get_radius_new(psf_pow, 2)

# show PSF
# img = Image_Plot()
# img.subplots(1,2)
# img.axs[0][0].set_title("PSF")
# fig = img.axs[0][0].imshow(psf_img)
# img.figure.colorbar(fig, ax=img.axs[0][0])
# img.axs[0][1].set_title("e-PSF")
# fig = img.axs[0][1].imshow(epsf_img)
# img.figure.colorbar(fig, ax=img.axs[0][1])
# img.show_img()

task_list = [[i,j] for i in range(si_num) for j in range(sr_num)]
my_task_list = tool_box.alloc(task_list, numprocs)[rank]

# simulate galaxy images

for ij in my_task_list:
    i,j = ij
    # measured shear
    mgs = numpy.zeros((2, shear_num))

    for k in range(shear_num):

        gal = galsim.Sersic(scale_radius=scale_radius[j], n=sersic_index[i], trunc=4.5 * scale_radius[j], flux=1.0)

        GN = numpy.zeros((3,))

        t1 = time.time()
        for ik in range(gal_num):

            gal_e = gal.shear(e1=gal_e1[ik], e2=gal_e2[ik])

            for ikk in range(4):
                #                     print(i,j,k,ik,ikk)
                gal_r = gal_e.rotate(theta=ikk * 45 * galsim.degrees)
                gal_s = gal_r.shear(g1=g1[k], g2=g2[k])
                # circle PSF
                gal_c = galsim.Convolve([gal_s, psf])
                gal_img = galsim.ImageD(stamp_size, stamp_size)
                gal_c.drawImage(image=gal_img, scale=pixel_scale)
                gal_img_arr = gal_img.array
                gal_pow = fq.pow_spec(gal_img_arr)

                mg1, mg2, mn = fq.shear_est(gal_pow, psf_pow, F=True)[:3]
                GN[0] += mg1
                GN[1] += mg2
                GN[2] += mn

        mgs[0, k] = GN[0] / GN[2]
        mgs[1, k] = GN[1] / GN[2]
        t2 = time.time()
        print(t2-t1)
    # if rank == 0:
    #     print(g1-mgs[0])
    #     print(g2-mgs[1])

    mc1 = tool_box.fit_1d(g1, mgs[0], 1, method="scipy")
    mc2 = tool_box.fit_1d(g2, mgs[1], 1, method="scipy")

    mcs[0, i, j] = mc1[1] -1 # m1
    mcs[1, i, j] = mc1[0] # c1
    mcs[2, i, j] = mc2[1]-1 # m2
    mcs[3, i, j] = mc2[0] # c2


if rank == 0:
    numpy.savez("./cache_%s.npz"%psf_tag, mcs)

    y = numpy.zeros((si_num*sr_num,))
    x = numpy.zeros((si_num*sr_num,))
    for i in range(si_num):
        for j in range(sr_num):
            tag = int(i*sr_num + j)
            y[tag] = sersic_index[i]
            x[tag] = scale_radius[j]

    titles = ["$m_1$","$c_1$","$m_2$", "$c_2$"]
    img = Image_Plot(xpad=0.25,ypad=0.25,fig_x=4,fig_y=3)
    img.subplots(2, 2)
    for i in range(2):
        for j in range(2):
            tag = int(i*2 + j)
            z = mcs[tag].flatten()
            img.scatter_pts(i,j,x,y,z,pts_size=30)

            img.set_label(i,j,0,"Sersic Index")
            img.set_label(i,j,1,"Scale Radius")
            img.axs[i][j].set_title(titles[tag])
    img.save_img("./mc_%s.png"%psf_tag)
    # img.show_img()