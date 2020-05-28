import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
import matplotlib.pyplot as plt
from matplotlib import cm


psf_type = "Moffat"
psf_flux = 1
psf_scale = 4
stamp_size = 48
seed = 56525

fq = Fourier_Quad(stamp_size, seed)
img = Image_Plot()
img.subplots(1, 1)
fig = img.axs[0][0].imshow(fq.kx2+fq.ky2)
img.figure.colorbar(fig, ax=img.axs[0][0])
img.axs[0][0].set_title("$k^2$")
img.del_ticks(0, 0, [0, 1])
img.save_img("E:/kxy.png")
img.show_img()
pst_num = 50


gal_fluxs = [4000, 16000, 32000, 100000]
steps = [0.6, 1, 2, 4]

fq = Fourier_Quad(stamp_size, seed)

max_radius = 7
pts = fq.ran_pts(pst_num, max_radius, step=1)
print(pts)
gal_img = fq.convolve_psf(pts, psf_scale, gal_fluxs[2] / pst_num, psf_type)

noise_1 = fq.draw_noise(0, 1)
noise_2 = fq.draw_noise(0, 1)
pnoise1 = fq.pow_spec(noise_1)
pnoise1_sqrt = numpy.sqrt(pnoise1)
pnoise2 = fq.pow_spec(noise_2)
noise_pow_diff = pnoise1 - pnoise2

gal_pow = fq.pow_spec(gal_img)
gal_pow_sqrt = numpy.sqrt(gal_pow)
galn_img = gal_img + noise_1
galn_pow = fq.pow_spec(galn_img)
cross_term = galn_pow - gal_pow - pnoise1

phase_gal = fq.pow_arg(gal_img)
phase_noise = fq.pow_arg(noise_1)
cos_2theta = 2*numpy.cos(phase_gal-phase_noise)

cross_term_rec = gal_pow_sqrt*pnoise1_sqrt*cos_2theta
plt_data = [[gal_pow_sqrt, pnoise1_sqrt, cos_2theta], [cross_term_rec,cross_term, cross_term - cross_term_rec]]
# plt_data = [numpy.sqrt(gal_pow)*numpy.sqrt(pnoise1)*2*numpy.cos(phase_gal - phase_noise),cross_term,]

img = Image_Plot()
img.subplots(2, 3)
for i in range(2):
    for j in range(3):
        fig = img.axs[i][j].imshow(plt_data[i][j])
        img.figure.colorbar(fig, ax=img.axs[i][j])
        img.del_ticks(i, j, [0, 1])
img.show_img()

exit()


img = Image_Plot()
img.subplots(4, 4)

theta = numpy.linspace(0, numpy.pi*2, 1000)
cos_theta = numpy.cos(theta)
sin_theta = numpy.sin(theta)

pkt_3d_data = [[],[]]
pow_angles = [[], [], []]

for k in range(4):
    gal_flux = gal_fluxs[3]

    fq = Fourier_Quad(stamp_size, seed)

    max_radius = 7
    pts = fq.ran_pts(pst_num, max_radius, step=steps[k])
    print(pts)
    gal_img = fq.convolve_psf(pts, psf_scale, gal_flux/pst_num, psf_type)

    noise_1 = fq.draw_noise(0, 1)
    noise_2 = fq.draw_noise(0, 1)
    pnoise1 = fq.pow_spec(noise_1)
    pnoise2 = fq.pow_spec(noise_2)
    noise_pow_diff = pnoise1 - pnoise2

    gal_pow = fq.pow_spec(gal_img)
    galn_img = gal_img + noise_1
    galn_pow = fq.pow_spec(galn_img)
    cross_term = galn_pow - gal_pow - pnoise1

    plt_data = [gal_img, cross_term, noise_pow_diff*numpy.sqrt(gal_pow)]

    phase_gal = fq.pow_arg(gal_img)
    phase_noise = fq.pow_arg(noise_1)
    pow_angles[0].append(phase_gal)
    pow_angles[1].append(phase_noise)
    pow_angles[2].append(phase_gal - phase_noise)

    pkt_3d_data[0].append(cross_term)
    pkt_3d_data[1].append(noise_pow_diff*numpy.sqrt(gal_pow))

    for j in range(3):
        fig = img.axs[k][j].imshow(plt_data[j])
        img.figure.colorbar(fig, ax=img.axs[k][j])
        img.del_ticks(k,j,[0,1])
    img.axs[k][3].scatter(pts[0], pts[1])
    img.axs[k][3].set_ylim(-max_radius*1.2,max_radius*1.2)
    img.axs[k][3].set_xlim(-max_radius*1.2,max_radius*1.2)
    img.axs[k][3].plot(max_radius*cos_theta, max_radius*sin_theta)

img.save_img("E:/gal_size.png")
img.show_img()
img.close_img()



for i in range(2):
    fig = plt.figure(figsize=(24, 5))
    axs = [fig.add_subplot(1, 4, i + 1 , projection='3d') for i in range(4)]
    for k in range(4):
        surf = axs[k].plot_surface(fq.mx, fq.my, pkt_3d_data[i][k], rstride=1, cstride=1, cmap=cm.coolwarm,
                                   linewidth=0, antialiased=False)
        fig.colorbar(surf, ax=axs[k], shrink=0.5, aspect=10)

    plt.savefig("E:/cross_term_%d.png"%i)
    # plt.show()
    plt.close()

img = Image_Plot()
img.subplots(3, 4)
for i in range(3):
    for j in range(4):
        fig = img.axs[i][j].imshow(pow_angles[i][j])
        img.figure.colorbar(fig, ax=img.axs[i][j])
        img.del_ticks(i,j,[0,1])
img.save_img("E:/arg.png")
