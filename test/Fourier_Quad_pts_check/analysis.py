import numpy
from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import h5py
from plot_tool import Image_Plot
import tool_box



def min_2d(coeff):
    # fxy = a + bx + cy + dx^2 + exy + fy^2
    a,b,c,d,e,f = coeff[0],coeff[1],coeff[2],coeff[3],coeff[4],coeff[5]
    det_a = 4*d*f-e*e
    x1_min = (-2*f*b+e*c)/det_a
    x1_sig = 1./numpy.sqrt(2*d)
    x2_min = (b*e-2*d*c)/det_a
    x2_sig = 1. / numpy.sqrt(2*f)
    return x1_min, x1_sig, x2_min, x2_sig



total_path = "E:/data/new_pdf/epsf/result_stack/weighted_u"

h5f = h5py.File(total_path + "/shear.hdf5", "r")
g1t = h5f["/g1"][()]
g2t = h5f["/g2"][()]
h5f.close()
print(g1t,g2t)
shear_num = g1t.shape[0]
g1m = numpy.zeros((4,shear_num))
g2m = numpy.zeros((4,shear_num))


chi_thresh = 200

h5f = h5py.File("%s/new_pdf_0.0030_0.0070_noise_free.hdf5"%total_path,"r")

shear_m = h5f["/PDF_new/shear"][()]

pic1_nm = total_path + "/result_0.0030_0.0070_noise_free_%d.png"%chi_thresh
pic2_nm = total_path + "/total_chisq_0.0030_0.0070_noise_free_%d_chisq.png"%chi_thresh
check_tag = 5

check = h5f["/PDF_new/%d/g1/check"%check_tag][()]

grid_row, grid_col = h5f["/PDF_new/%d/g1/chisq_grid"%check_tag][()].shape

total_chisq1 = numpy.zeros((grid_row*2,grid_col*5))
total_chisq2 = numpy.zeros((grid_row*2,grid_col*5))

for tag in range(10):
    m, n = divmod(tag, 5)

    g1_grid1 = h5f["/PDF_new/%d/g1/grid1"%tag][()]
    g1_grid2 = h5f["/PDF_new/%d/g1/grid2"%tag][()]
    chisq1 = h5f["/PDF_new/%d/g1/chisq_grid"%tag][()]

    g2_grid1 = h5f["/PDF_new/%d/g2/grid1"%tag][()]
    g2_grid2 = h5f["/PDF_new/%d/g2/grid2"%tag][()]
    chisq2 = h5f["/PDF_new/%d/g2/chisq_grid"%tag][()]

    check1_pdf_ori = h5f["/PDF_new/%d/g1/check"%tag][()]
    check2_pdf_ori = h5f["/PDF_new/%d/g2/check"%tag][()]

    total_chisq1[grid_row*m:grid_row*(m+1),grid_col*n:grid_col*(n+1)] = chisq1
    total_chisq2[grid_row*m:grid_row*(m+1),grid_col*n:grid_col*(n+1)] = chisq2

    idx1 = chisq1 < chi_thresh
    coeff_py1 = tool_box.fit_2d(g1_grid1[idx1], g1_grid2[idx1], chisq1[idx1],2)[0]
    x1,x1_sig, x2, x2_sig = min_2d(coeff_py1)
    g1m[0,tag] = x1
    g1m[1,tag] = x1_sig
    g1m[2,tag] = x2
    g1m[3,tag] = x2_sig

    idx2 = chisq2 < chi_thresh
    coeff_py2 = tool_box.fit_2d(g2_grid1[idx2], g2_grid2[idx2], chisq2[idx2],2)[0]
    x1, x1_sig, x2, x2_sig = min_2d(coeff_py2)
    g2m[0,tag] = x1
    g2m[1,tag] = x1_sig
    g2m[2,tag] = x2
    g2m[3,tag] = x2_sig

    text_str = "true g1: %.5f  g2: %.5f \n" \
               "g1_N:  %.5f(%.5f) g2_N: %.5f(%.5f)  " \
               "g1_U:  %.5f(%.5f) g2_U: %.5f(%.5f)"%\
               (g1t[tag],g2t[tag],g1m[0,tag],g1m[1,tag],g2m[0,tag],g2m[1,tag],
                g1m[2,tag],g1m[3,tag],g2m[2,tag],g2m[3,tag])
    #

    img = Image_Plot(fig_x=3, fig_y=3, ypad=0.15, xpad=0.1)
    img.subplots(2, 3)
    fig = img.axs[0][0].imshow(g1_grid1)
    img.figure.colorbar(fig, ax=img.axs[0][0],orientation='horizontal')
    img.axs[0][0].set_title("$g1_N grid$")

    fig = img.axs[0][1].imshow(g1_grid2)
    img.figure.colorbar(fig, ax=img.axs[0][1],orientation='horizontal')
    img.axs[0][1].set_title("$g1_U grid$")

    fig = img.axs[0][2].imshow(chisq1)
    img.figure.colorbar(fig, ax=img.axs[0][2],orientation='horizontal')
    img.axs[0][2].set_title("$g1 \chi^2 grid$")

    fig = img.axs[1][0].imshow(g2_grid1)
    img.figure.colorbar(fig, ax=img.axs[1][0],orientation='horizontal')
    img.axs[1][0].set_title("$g2_N grid$")

    fig = img.axs[1][1].imshow(g2_grid2)
    img.figure.colorbar(fig, ax=img.axs[1][1],orientation='horizontal')
    img.axs[1][1].set_title("$g2_U grid$")

    fig = img.axs[1][2].imshow(chisq2)
    img.figure.colorbar(fig, ax=img.axs[1][2],orientation='horizontal')
    img.axs[1][2].set_title("$g2 \chi^2 grid$")

    img.axs_text(0, 0, 1.2, 0.5, text_str,text_fontsize=10)
    img.save_img(total_path + "/shear_point_chisq_0.0030_0.0070_noise_free_%d.png" % tag)
    # img.show_img()
    img.close_img()



mc1 = numpy.array(tool_box.data_fit(g1t, g1m[0], g1m[1]))
mc2 = numpy.array(tool_box.data_fit(g2t, g2m[0], g2m[1]))
mc1[0] = mc1[0] - 1
mc2[0] = mc2[0] - 1
print(mc1)
print(mc2)
text_str = "g_N\nm1: %.5f(%.5f)\nc1: %.5f(%.5f)\n" \
           "m2: %.5f(%.5f)\nc2: %.5f(%.5f)\n"\
           %(mc1[0],mc1[1],mc1[2],mc1[3],mc2[0],mc2[1],mc2[2],mc2[3])

mc1_u = numpy.array(tool_box.data_fit(g1t, g1m[2], g1m[3]))
mc2_u = numpy.array(tool_box.data_fit(g2t, g2m[2], g2m[3]))
mc1_u[0] = mc1_u[0] - 1
mc2_u[0] = mc2_u[0] - 1
print(mc1_u)
print(mc2_u)
#
text_str += "\ng_U\nm1: %.5f(%.5f)\nc1: %.5f(%.5f)\n" \
           "m2: %.5f(%.5f)\nc2: %.5f(%.5f)\n"\
            %(mc1_u[0],mc1_u[1],mc1_u[2],mc1_u[3],mc2_u[0],mc2_u[1],mc2_u[2],mc2_u[3])


mc1 = numpy.array(tool_box.data_fit(g1t, shear_m[0], shear_m[1]))
mc2 = numpy.array(tool_box.data_fit(g2t, shear_m[2], shear_m[3]))
mc1[0] = mc1[0] - 1
mc2[0] = mc2[0] - 1
text_str += "\nori_PDF\nm1: %.5f(%.5f)\nc1: %.5f(%.5f)\n" \
           "m2: %.5f(%.5f)\nc2: %.5f(%.5f)\n"\
            %(mc1[0],mc1[1],mc1[2],mc1[3],mc2[0],mc2[1],mc2[2],mc2[3])

mc1 = numpy.array(tool_box.data_fit(g1t, shear_m[4], shear_m[5]))
mc2 = numpy.array(tool_box.data_fit(g2t, shear_m[6], shear_m[7]))
mc1[0] = mc1[0] - 1
mc2[0] = mc2[0] - 1
text_str += "\nAverage\nm1: %.5f(%.5f)\nc1: %.5f(%.5f)\n" \
           "m2: %.5f(%.5f)\nc2: %.5f(%.5f)\n"\
            %(mc1[0],mc1[1],mc1[2],mc1[3],mc2[0],mc2[1],mc2[2],mc2[3])


img = Image_Plot(pts_size=3,fig_x=8,fig_y=6)
img.subplots(1,1)
img.axs[0][0].plot([-0.045,0.045],[-0.045,0.045], c="grey",alpha=0.3, ls="dotted",label="y=x")
img.axs[0][0].errorbar(g1t, g1m[0], g1m[1], fmt=" ", label="g1_N",capsize=2, marker="s",ms=img.pts_size)
img.axs[0][0].errorbar(g2t, g2m[0], g2m[1], fmt=" ", label="g2_N",capsize=2, marker="s",ms=img.pts_size)
img.axs[0][0].errorbar(g1t, g1m[2], g1m[3], fmt=" ", label="g1_U",capsize=2, marker="s",ms=img.pts_size)
img.axs[0][0].errorbar(g2t, g2m[2], g2m[3], fmt=" ", label="g2_U",capsize=2, marker="s",ms=img.pts_size)
img.axs[0][0].errorbar(g1t, shear_m[0], shear_m[1], fmt=" ", label="g1_pdf_ori",capsize=2, marker="s",ms=img.pts_size)
img.axs[0][0].errorbar(g2t, shear_m[2], shear_m[3], fmt=" ", label="g2_pdf_ori",capsize=2, marker="s",ms=img.pts_size)
img.axs_text(0,0, 0.55, 0.01,text_str,text_fontsize=12)
img.axs[0][0].legend(loc="lower right")
img.axs[0][0].set_xlim(-0.045,0.045)
img.axs[0][0].set_ylim(-0.045,0.045)
img.save_img(pic1_nm)
# img.show_img()


idx1 = total_chisq1 > chi_thresh
total_chisq1[idx1] = numpy.nan
idx2 = total_chisq2 > chi_thresh
total_chisq2[idx2] = numpy.nan

img = Image_Plot(fig_x=10, fig_y=4)
img.subplots(2, 1)
fig = img.axs[0][0].imshow(total_chisq1)
img.figure.colorbar(fig, ax=img.axs[0][0])

fig = img.axs[1][0].imshow(total_chisq2)
img.figure.colorbar(fig, ax=img.axs[1][0])
img.save_img(pic2_nm)
# img.show_img()

