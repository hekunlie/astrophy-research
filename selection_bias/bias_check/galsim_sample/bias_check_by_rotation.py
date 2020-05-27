import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
from sys import path, argv
path.append('%s/work/mylib/' % my_home)
import numpy
import galsim
from astropy.io import fits
from Fourier_Quad import Fourier_Quad
from mpi4py import MPI
import h5py
import tool_box
from plot_tool import Image_Plot


#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
numprocs = comm.Get_size()


grid_num = 30

itemsize = MPI.DOUBLE.Get_size()
element_num = grid_num*grid_num*4
if rank == 0:
    # bytes for 10 double elements
    nbytes = element_num*itemsize
else:
    nbytes = 0

# on rank 0 of comm, create the contiguous shared block
win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
result_buff = numpy.ndarray(buffer=buf1, dtype='d', shape=(grid_num,grid_num*4)) # array filled with zero

total_grid_num = grid_num*grid_num

my_grids = tool_box.alloc([i for i in range(total_grid_num)], numprocs)[rank]

idxs = numpy.linspace(0.31, 1.2, grid_num)
radii = numpy.linspace(0.2, 0.8, grid_num)

gal_flux = 6000
noise_sig = 60
shear_num = 40

pixel_scale = 0.187
stamp_size = 64

g1_result = numpy.zeros((2, shear_num))
g2_result = numpy.zeros((2, shear_num))

temp = numpy.zeros((4, 3))

fq = Fourier_Quad(stamp_size, 123)

psf = galsim.Moffat(beta=3.5, fwhm=0.7, flux=1.0, trunc=1.4).shear(e1=0.15, e2=0.)

psf_img = galsim.ImageD(stamp_size, stamp_size)
psf.drawImage(image=psf_img, scale=pixel_scale)
psf_arr = psf_img.array
# hdu = fits.PrimaryHDU(psf_arr)
# psf_path = './psf.fits'
# hdu.writeto(psf_path, overwrite=True)

psf_pow = fq.pow_spec(psf_arr)
fq.get_radius_new(psf_pow, 2)

shear_path = "./shear.hdf5"
h5f = h5py.File(shear_path, "r")
g1_input = h5f["/g1"][()]
g2_input = h5f["/g2"][()]
h5f.close()


for ij in my_grids:
    m, n = divmod(ij, grid_num)
    sersic_idx = idxs[m]
    scale_radius = radii[n]

    # img_buf = numpy.zeros((4*stamp_size, shear_num*stamp_size))

    gal = galsim.Sersic(scale_radius=scale_radius, n=sersic_idx, trunc=4.5 * scale_radius, flux=1.0)

    e = numpy.random.uniform(0, 0.7, shear_num)
    theta = numpy.random.normal(0, numpy.pi * 2, shear_num)
    cos_2theta = numpy.cos(2 * theta)
    sin_2theta = numpy.sin(2 * theta)

    for i in range(shear_num):

        g1 = g1_input[i]
        g2 = g2_input[i]

        for j in range(4):

            e1 = e[i]*cos_2theta[i]
            e2 = e[i]*sin_2theta[i]
            # print(e1, e2)
            gal_e = gal.shear(e1=e1, e2=e2).withFlux(gal_flux).rotate(j*180./4*galsim.degrees)

            gal_s = gal_e.shear(g1=g1, g2=g2)
            gal_c = galsim.Convolve([gal_s, psf])
            gal_img = galsim.ImageD(stamp_size, stamp_size)
            gal_c.drawImage(image=gal_img, scale=pixel_scale)

            gal_arr = gal_img.array
            gal_pow = fq.pow_spec(gal_arr)

            temp[j] = fq.shear_est(gal_pow, psf_pow, F=True)[:3]

            # img_buf[j*stamp_size:(j+1)*stamp_size, i*stamp_size:(i+1)*stamp_size] = gal_img.array

        gh1, gh1_sig = fq.find_shear_mean(temp[:,0], temp[:,2])
        gh2, gh2_sig = fq.find_shear_mean(temp[:,1], temp[:,2])

        # print(g1, gh1,gh1_sig)
        # print(g2, gh2,gh2_sig)

        g1_result[0,i] = gh1
        g1_result[1,i] = gh1_sig

        g2_result[0,i] = gh2
        g2_result[1,i] = gh2_sig

    # hdu = fits.PrimaryHDU(img_buf)
    # psf_path = './gal_imgs.fits'
    # hdu.writeto(psf_path, overwrite=True)

    mc1 = numpy.array(tool_box.data_fit(g1_input, g1_result[0],g1_result[0]))
    mc2 = numpy.array(tool_box.data_fit(g2_input, g2_result[0],g2_result[0]))
    print(mc1)
    print(mc2)

    mc1[0] = mc1[0] - 1
    mc2[0] = mc2[0] - 1

    text_str = "m1: %.5f(%.5f) c1: %.5f(%.5f)\n"%(mc1[0], mc1[1], mc1[2], mc1[3]) + \
               "m2: %.5f(%.5f) c2: %.5f(%.5f)"%(mc2[0], mc2[1], mc2[2], mc2[3])
    print(text_str)

    result_buff[m, n] = mc1[0]
    result_buff[m, n + grid_num] = mc1[2]
    result_buff[m, n + grid_num*2] = mc2[0]
    result_buff[m, n + grid_num*3] = mc2[2]
    # img = Image_Plot(fig_x=6, fig_y=4,xpad=0.25)
    # img.subplots(1,2)
    # img.axs[0][0].errorbar(g1_input, g1_result[0], g1_result[1], c="C1", label="g1",fmt=" ")
    # img.axs[0][0].scatter(g1_input, g1_result[0],c="C1")
    # img.axs[0][0].plot(g1_input, g1_input,ls="--", c="grey")
    #
    # img.axs[0][0].errorbar(g2_input, g2_result[0], g2_result[1], c="C2", label="g2", fmt=" ")
    # img.axs[0][0].scatter(g2_input, g2_result[0], c="C2")
    #
    # img.axs[0][0].set_xlim(-0.06, 0.06)
    # img.axs[0][0].set_ylim(-0.06, 0.06)
    # img.axs[0][0].legend()
    #
    #
    # diff1 = g1_result[0] - g1_input
    # diff2 = g2_result[0] - g2_input
    # y1 = max([diff1.max(), diff2.max()])
    # y2 = min([diff1.min(), diff2.min()])
    # dy = (y1 - y2)*0.1
    #
    # img.axs[0][1].scatter(g1_input,diff1, label="$\delta$g1")
    # img.axs[0][1].scatter(g2_input,diff2, label="$\delta$g2")
    # img.axs[0][1].plot([-0.04,0.04], [0,0],ls="--", c="grey")
    # img.axs[0][1].set_ylim(y2-dy, y1+dy)
    #
    # img.axs[0][1].legend()
    #
    # img.axs_text(0,0, 0.9, 0.07, text_str, text_fontsize=12)
    # img.save_img("./img.png")
comm.Barrier()
if rank == 0:
    numpy.savez("./grid_data.npz", idxs, radii, result_buff)
comm.Barrier()