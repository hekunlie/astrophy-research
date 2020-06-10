from sys import path
path.append("D:/GitHub/astrophy-research/mylib")
import numpy
from Fourier_Quad import Fourier_Quad
import tool_box
from plot_tool import Image_Plot
import h5py

psf_type = "Moffat"
psf_flux = 1
psf_scale = 4
stamp_size = 44
seed = 2301652

pts_num = 40
max_radius = 7
gal_flux = 8000/pts_num


shear_path = "./shear.hdf5"
h5f = h5py.File(shear_path, "r")
g1_input = h5f["/g1"][()]
g2_input = h5f["/g2"][()]
h5f.close()

shear_num = g1_input.shape[0]

fq = Fourier_Quad(stamp_size, seed)

psf_img = fq.cre_psf(psf_scale, psf_flux, psf_type)
psf_pow = fq.pow_spec(psf_img)
fq.get_radius(psf_pow, 2)
print(fq.hlr)
img = Image_Plot()
img.subplots(1,2)
img.axs[0][0].imshow(psf_img)
img.axs[0][1].imshow(psf_pow)
img.show_img()

num_scale = 10
rotation = 4

temp = numpy.zeros((rotation*num_scale,3))

g1_result = numpy.zeros((2,shear_num))
g2_result = numpy.zeros((2,shear_num))

for i in range(shear_num):
    g1 = g1_input[i]
    g2 = g2_input[i]

    for j in range(num_scale):

        pts = fq.ran_pts(pts_num, max_radius)
        # print(pts)
        tag = j*rotation
        for k in range(rotation):

            pts_r = fq.rotate(pts, numpy.pi/rotation*k)

            pst_s = fq.shear(pts_r, g1, g2)

            noise_1 = fq.draw_noise(0, 60)
            # noise_2 = fq.draw_noise(0, 0.01)
            # pnoise = fq.pow_spec(noise_2)

            gal_img = fq.convolve_psf(pst_s, psf_scale, gal_flux, psf_type)+noise_1
            img = Image_Plot()
            img.subplots(1,2)
            idx = gal_img > 90
            gal_temp = numpy.zeros_like(gal_img)
            gal_temp[idx] = gal_img[idx]
            img.axs[0][0].imshow(gal_img)
            img.axs[0][1].imshow(gal_temp)
            img.show_img()
            gal_pow = fq.pow_spec(gal_img)

            temp[tag + k] = fq.shear_est(gal_pow, psf_pow, F=True)[:3]

    gh1, gh1_sig = fq.find_shear_mean(temp[:,0], temp[:,2])
    gh2, gh2_sig = fq.find_shear_mean(temp[:,1], temp[:,2])

    g1_result[0,i] = gh1
    g1_result[1,i] = gh1_sig

    g2_result[0,i] = gh2
    g2_result[1,i] = gh2_sig

mc1 = numpy.array(tool_box.data_fit(g1_input, g1_result[0],g1_result[0]))
mc2 = numpy.array(tool_box.data_fit(g2_input, g2_result[0],g2_result[0]))
print(mc1)
print(mc2)

mc1[0] = mc1[0] - 1
mc2[0] = mc2[0] - 1

text_str = "m1: %.5f(%.5f) c1: %.5f(%.5f)\n"%(mc1[0], mc1[1], mc1[2], mc1[3]) + \
           "m2: %.5f(%.5f) c2: %.5f(%.5f)"%(mc2[0], mc2[1], mc2[2], mc2[3])
print(text_str)

img = Image_Plot(fig_x=6, fig_y=4,xpad=0.25)
img.subplots(1,2)
img.axs[0][0].errorbar(g1_input, g1_result[0], g1_result[1], c="C1", label="g1",fmt=" ")
img.axs[0][0].scatter(g1_input, g1_result[0],c="C1")
img.axs[0][0].plot(g1_input, g1_input,ls="--", c="grey")

img.axs[0][0].errorbar(g2_input, g2_result[0], g2_result[1], c="C2", label="g2", fmt=" ")
img.axs[0][0].scatter(g2_input, g2_result[0], c="C2")

img.axs[0][0].set_xlim(-0.06, 0.06)
img.axs[0][0].set_ylim(-0.06, 0.06)
img.axs[0][0].legend()
img.axs_text(0,0, 0.9, 0.07, text_str, text_fontsize=12)


diff1 = g1_result[0] - g1_input
diff2 = g2_result[0] - g2_input
y1 = max([diff1.max(), diff2.max()])
y2 = min([diff1.min(), diff2.min()])
dy = (y1 - y2)*0.1

img.axs[0][1].scatter(g1_input,diff1, label="$\delta$g1")
img.axs[0][1].scatter(g2_input,diff2, label="$\delta$g2")
img.axs[0][1].plot([-0.04,0.04], [0,0],ls="--", c="grey")
img.axs[0][1].set_ylim(y2-dy, y1+dy)

img.axs[0][1].legend()
img.save_img("./img.png")
img.show_img()
