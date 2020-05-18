from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot


psf_type = "GAUSS"
psf_flux = 1
psf_scale = 4
stamp_size = 64
seed = 2301

pst_num = 100
max_radius = 8
gal_flux = 100


g1, g2 = 0.03555, -0.0433

fq = Fourier_Quad(stamp_size, seed)

psf_img = fq.cre_psf(psf_scale, psf_flux, psf_type)
psf_pow = fq.pow_spec(psf_img)
fq.get_radius_new(psf_pow, 2)
print(fq.hlr)

pts = fq.ran_pts(pst_num, max_radius)

num_scale = 1
rotation = 4
mg1 = numpy.zeros((rotation*num_scale,))
mg2 = numpy.zeros((rotation*num_scale,))
mn = numpy.zeros((rotation*num_scale,))
mu = numpy.zeros((rotation*num_scale,))

img = Image_Plot()
img.subplots(2,5)
img.axs[0][0].imshow(psf_img)

for j in range(num_scale):
    tag = j*rotation
    for i in range(rotation):

        pts_r = fq.rotate(pts, numpy.pi/rotation*i)

        pst_s = fq.shear(pts_r, g1, g2)

        noise_1 = fq.draw_noise(0, 0.01)
        noise_2 = fq.draw_noise(0, 0.01)
        pnoise = fq.pow_spec(noise_2)

        gal_img = fq.convolve_psf(pst_s, psf_scale, gal_flux, psf_type)
        gal_pow = fq.pow_spec(gal_img)

        mg1[tag+i],mg2[tag+i],mn[tag+i],mu[tag+i] = fq.shear_est(gal_pow, psf_pow, F=True)[:4]

        # img.axs[0][i+1].imshow(gal_img)
        # img.axs[1][i+1].imshow(gal_img/peak)



est_g1, est_g2 = mg1.sum()/mn.sum(), mg2.sum()/mn.sum()
# print(g1/est_g1,g2/est_g2/numpy.sqrt(2))
print("The est  g1: %10.5f, g2: %10.5f"%(est_g1, est_g2))
print("The true g1: %10.5f, g2: %10.5f"%(g1, g2))
print("The diff g1: %10.5f, g2: %10.5f"%(g1-est_g1, g2-est_g2))
# img.show_img()
img.close_img()
print(mg1)
print(mg2)
print(mn)