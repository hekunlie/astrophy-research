from Fourier_Quad import *
import numpy
import matplotlib.pyplot as plt
import galsim

og1 = numpy.linspace(-0.01,0.01,11)
og2 = numpy.linspace(-0.01,0.01,11)
numpy.random.shuffle(og2)

size = 60
pixel_scale = 0.3
psf_ori = galsim.Exponential(half_light_radius=1.2)
psf_img = psf_ori.drawImage(nx=size,ny=size,scale=pixel_scale)
psf_arr = psf_img.array
gal_num =200000
for i in range(11):
    g1 = og1[i]
    g2 = og2[i]
    data_matrix = 0
    p=0
    while True:
        oe = numpy.random.normal(loc=0,scale=0.15,size=gal_num)
        theta = 4*numpy.random.uniform(0,1,size=gal_num)*numpy.pi
        oe1 = oe*numpy.cos(theta)
        oe2 = oe*numpy.sin(theta)
        if numpy.abs(numpy.mean(oe1))<1e-5 and numpy.abs(numpy.mean(oe2))<1e-5:
            break
    print(numpy.mean(oe1),numpy.mean(oe2))

    for m in range(gal_num):
        e1 = oe1[m]
        e2 = oe2[m]
        gal_o = galsim.Gaussian(half_light_radius=2.2,flux = 1.e6)
        gal_e = gal_o.shear(e1=e1,e2=e2)
        gal_s = gal_e.shear(g1=g1,g2=g2)
        gal_c = galsim.Convolve([gal_s,psf_ori])
        final = gal_c.drawImage(nx=size,ny=size,scale=pixel_scale).array
        G1,G2,N,U = Fourier_Quad().shear_est(final,psf_arr,size,F=False)[0:4]
        ith_row = numpy.array([ G1, G2, N,U,g1,g2])
        if p == 0:
            data_matrix = ith_row
            p=1
        else:
            data_matrix = numpy.row_stack((data_matrix, ith_row))
    cache_name = 'gal_test_%.3f_%.3f'%(g1,g2)
    numpy.savez(cache_name,data_matrix)
    G1 = data_matrix[:,0]
    G2 = data_matrix[:,1]
    N = data_matrix[:,2]
    U = data_matrix[:,3]
    g1 = data_matrix[0,-2]
    g2 = data_matrix[0,-1]
    g1_est_mean = numpy.sum(G1)/numpy.sum(N)
    g1_est_mean_sig = numpy.std(G1/N)/numpy.sqrt(len(N))
    g2_est_mean = numpy.sum(G2) / numpy.sum(N)
    print(g1,g1_est_mean,g2,g2_est_mean)