import matplotlib
matplotlib.use("Agg")
import numpy
# import os
# my_home = os.popen("echo $HOME").readlines()[0][:-1]
# from sys import path
# path.append('%s/work/fourier_quad/'%my_home)
# # path.append("E:/Github/astrophy-research/")
# import time
# from Fourier_Quad import Fourier_Quad
# # import galsim
# import matplotlib.pyplot as plt
# from astropy.io import fits
# import tool_box
# from mpi4py import MPI
# #
#

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()

# if rank < 10:
#     g = numpy.load("./cor_g.npz")['arr_0'][:, 0]
# else:
#     g = numpy.load("./cor_g.npz")['arr_0'][:, 1]

ig1 = numpy.load("./g.npz")['arr_0'][:,0]
ig2 = numpy.load("./g.npz")['arr_0'][:,1]

num = 200000
stamp_size = 54
pixel_scale = 0.2
fq = Fourier_Quad(stamp_size, rank*10+111)
psf_img = fq.cre_psf(4,1,"Moffat")
psf_ps = fq.pow_spec(psf_img)

# psf = galsim.Moffat(beta=3.5, scale_radius=1.0, flux=1.0, trunc=3)
# psf_img = galsim.ImageD(stamp_size, stamp_size)
# psf.drawImage(image=psf_img, scale=pixel_scale)
# psf_img = psf_img.array

data = numpy.zeros((num,4))
gal_pool = []
noise_sig = 1
ra = numpy.random.uniform(0.6, 1.2, num)
g1,g2 = ig1[rank], ig2[rank]
t1 = time.time()
psf_res = fq.get_radius_new(psf_ps, 2)
for i in range(num):
    pts = fq.ran_pos(45,8,(ig1[rank],ig2[rank]))[1]
    gal_img = fq.convolve_psf(pts,4,300,"Moffat") + fq.draw_noise(0,noise_sig)
    noise = fq.draw_noise(0,noise_sig)
    #
    # gal = galsim.Sersic(half_light_radius=ra[i], n=3, trunc=4.5*ra[i])
    # gal = galsim.Gaussian(flux=1.e4, half_light_radius=ra[i])
    # for j in range(4):
    #     gal_s = gal.rotate(j/4.*numpy.pi*galsim.radians)
    #     gal_g = gal_s.shear(g1=g1, g2=g2)
    #     gal_c = galsim.Convolve([gal_g, psf])
    #     img = galsim.ImageD(stamp_size, stamp_size)
    #     gal_c.drawImage(image=img, scale=pixel_scale)
    #     gal_img = img.array# + fq.draw_noise(0, noise_sig)
    #     # noise = fq.draw_noise(0,noise_sig)

    res = fq.shear_est(gal_img, psf_ps, noise, F=True)
    if rank == 0 and i < 100:
        gal_pool.append(gal_img)
    data[i] = res[0],res[1],res[2],res[3]
t2 = time.time()
if rank == 0:
    print(t2-t1)
    img = fq.stack(gal_pool,10)
    img_w = fits.PrimaryHDU(img)
    img_w.writeto("./gals.fits",overwrite=True)
    numpy.savez("./data.npz",data)
MG1 = data[:, 0]
MG2 = data[:, 1]
MN = data[:, 2]
MU = data[:, 3]

# be careful that the "MU" defined in FRESH is the different from that in ours
# MN + (-) MU for our definition (g1(2)) of MU and MV which is the same as
# those in the paper Zhang et al. 2017 ApJ, 834:8
DE1 = MN + MU
DE2 = MN - MU

g1_h, g1_sig = fq.fmin_g(MG1, DE1, bin_num=8)
g2_h, g2_sig = fq.fmin_g(MG2, DE2, bin_num=8)

g1_ho, g1_sigo = fq.fmin_g_new(MG1, DE1, bin_num=8)
g2_ho, g2_sigo = fq.fmin_g_new(MG2, DE2, bin_num=8)

g1_hm = MG1.mean()/MN.mean()
g1_sigm = numpy.sqrt((MG1**2).mean()/MN.mean()**2)/numpy.sqrt(num)
g2_hm = MG2.mean()/MN.mean()
g2_sigm = numpy.sqrt((MG2**2).mean()/MN.mean()**2)/numpy.sqrt(num)

result = [g1_h, g1_sig, g2_h, g2_sig, g1_ho, g1_sigo, g2_ho, g2_sigo, g1_hm, g1_sigm, g2_hm, g2_sigm]
gs = comm.gather(result, root=0)

if rank == 0:

    gs = numpy.array(gs)

    g1 = gs[:,0]
    dg1 = gs[:,1]
    g2 = gs[:,2]
    dg2 = gs[:,3]
    emc1 = tool_box.data_fit(ig1, g1, dg1)
    emc2 = tool_box.data_fit(ig2, g2, dg2)

    g1o = gs[:,4]
    dg1o = gs[:,5]
    g2o = gs[:,6]
    dg2o = gs[:,7]
    emc1o = tool_box.data_fit(ig1, g1o, dg1o)
    emc2o = tool_box.data_fit(ig2, g2o, dg2o)

    g1m = gs[:,8]
    dg1m = gs[:,9]
    g2m = gs[:,10]
    dg2m = gs[:,11]
    emc1m = tool_box.data_fit(ig1, g1m, dg1m)
    emc2m = tool_box.data_fit(ig2, g2m, dg2m)

    plt.subplot(121)
    x = numpy.linspace(-0.03,0.03,5)
    plt.errorbar(ig1,g1,dg1,capsize=3,fmt="none")
    plt.plot(x, emc1[0]*x+emc1[2])
    plt.subplot(122)
    plt.errorbar(ig2,g2,dg2,capsize=3,fmt="none")
    plt.plot(x, emc2[0]*x+emc2[2])
    plt.savefig("./test.png")
    # plt.show()
    print("m: %8.5f (%10.7f) c: %8.5f (%10.7f)"%(emc1[0]-1,emc1[1],emc1[2],emc1[3]))
    print("m: %8.5f (%10.7f) c: %8.5f (%10.7f)"%(emc2[0]-1,emc2[1],emc2[2],emc2[3]))
    print("m: %8.5f (%10.7f) c: %8.5f (%10.7f)"%(emc1o[0]-1,emc1o[1],emc1o[2],emc1o[3]))
    print("m: %8.5f (%10.7f) c: %8.5f (%10.7f)"%(emc2o[0]-1,emc2o[1],emc2o[2],emc2o[3]))
    print("m: %8.5f (%10.7f) c: %8.5f (%10.7f)"%(emc1m[0]-1,emc1m[1],emc1m[2],emc1m[3]))
    print("m: %8.5f (%10.7f) c: %8.5f (%10.7f)"%(emc2m[0]-1,emc2m[1],emc2m[2],emc2m[3]))
#




