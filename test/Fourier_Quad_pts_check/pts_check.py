from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot


psf_type = "Moffat"
psf_flux = 1
psf_scale = 4
stamp_size = 64
seed = 12322

pst_num = 100
max_radius = 10
gal_flux = 1000

g1, g2 = -0.01, 0.03

fq = Fourier_Quad(stamp_size, seed)

psf_img = fq.cre_psf(psf_scale, psf_flux, psf_type)

fq.get_radius_new(psf_img,2)

img = Image_Plot()
img.subplots(1,5)
img.axs[0][0].imshow(psf_img)

pts = fq.ran_pts(pst_num, max_radius)

shear_estimator = numpy.zeros((4, 5))

for i in range(4):
    pts_r = fq.rotate(pts,numpy.pi/4*i)
    pst_s = fq.shear(pts_r, g1, g2)

    gal_img = fq.convolve_psf(pst_s, psf_scale, gal_flux, psf_type)

    shear_estimator[i] = fq.shear_est(gal_img, psf_img)

    img.axs[0][i+1].imshow(gal_img)

img.show_img()
img.close_img()

est_g1, est_g2 = shear_estimator[:,0].mean()/shear_estimator[:,2].mean(), shear_estimator[:,1].mean()/shear_estimator[:,2].mean()
print("The est  g1: %10.5f, g2: %10.5f"%(est_g1, est_g2))
print("The true g1: %10.5f, g2: %10.5f"%(g1, g2))
