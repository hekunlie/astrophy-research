import numpy
from numpy import fft
from scipy.optimize import least_squares
from sys import path
path.append("E:/Github/astrophy-research/mylib")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import tool_box
import h5py
from astropy.cosmology import FlatLambdaCDM
from astropy.cosmology import FlatwCDM
from mpl_toolkits.axes_grid1.inset_locator import zoomed_inset_axes
from mpl_toolkits.axes_grid1.inset_locator import mark_inset



def get_dist(H0, Om0, z, w=None):
    h = H0/100
    if w:
        cosmos = FlatwCDM(H0=H0, Om0=Om0, w0=w)
    else:
        cosmos = FlatLambdaCDM(H0=H0, Om0=Om0)
    return cosmos.comoving_distance(z).value*h

def gamma_coeff(z_lens, delta_sigma, H0, Om0, w=None):
    dist = get_dist(H0, Om0, z_lens, w)
    return delta_sigma*dist*(1+z_lens)**2*6.013625/10**7


w = -0.8
npara = 5
num = 100
z_lens = 0.6
z_max = 1.2
Om0 = 0.25
H0 = 70
delta_sigma = 15

# img = Image_Plot(fig_x=8,fig_y=6)
# img.subplots(1, 3)
diff_ratio = numpy.zeros((10,))
z_fore = numpy.linspace(0.1, 1, 10)
for m in range(10):
    max_ratio = 0
    ratio_tag = -1
    z = numpy.linspace(z_fore[m], z_max, num)
    fore_tag = 0
    for k in range(5, num-5):
        len_tag1 = k
        dt = 0
        dist_c = numpy.zeros((npara, num))

        ratio = numpy.zeros((npara, num-len_tag1-dt))

        tag = tool_box.find_near(z[len_tag1 + dt:], z_max - 0.02)

        for i in range(npara):
            Om0_ = Om0 + i * 0.02
            # distance
            dist_lens_c = get_dist(H0, Om0_, z_lens, w)

            for j in range(num):
                dist_c[i, j] = get_dist(H0, Om0_, z[j], w)

            # print(numpy.max(numpy.abs(dist_c[i]/H0*100/(1+z) - dist_p[i])))

            gamma_z1 = 1 - dist_c[i, fore_tag]/dist_c[i, len_tag1]
            gamma_z2 = 1 - dist_c[i, fore_tag]/dist_c[i, len_tag1+dt:]
            ratio[i] = gamma_z1/gamma_z2

            # img.axs[0][0].plot(z[len_tag1+dt:],ratio[i],c="C%d"%i,label="$\Omega_{m0}=%.2f$"%Om0_, linewidth=2.5)
            #
            # img.axs[0][1].plot(z[len_tag1+dt:][tag:],ratio[i,tag:],c="C%d"%i,label="$\Omega_{m0}=%.2f$"%Om0_, linewidth=2.5)

        if numpy.max((ratio[0] - ratio[-1])*100)>= max_ratio:
            ratio_tag = k
            max_ratio = numpy.max((ratio[0] - ratio[-1])*100)
    print(ratio_tag, max_ratio,z[ratio_tag])
    diff_ratio[m] = max_ratio
img = Image_Plot(fig_x=8,fig_y=6)
img.subplots(1, 1)
img.axs[0][0].plot(z_fore, diff_ratio)
img.set_label(0,0,0,"R")
img.set_label(0,0,1,"Foreground Z")
img.save_img("E:/works/Galaxy-Galaxy_lensing/gamma_ratio/max_gamma_ratio.png")
img.show_img()

#
# img.axs[0][2].plot(z[len_tag1+dt:],(ratio[0] - ratio[-1])*100, label="$10^2[R(\Omega_{m0}=%.2f) -R(\Omega_{m0}=%.2f)]$"%(Om0, Om0 + (npara-1)* 0.02), linewidth=2.5)
#
# img.axs[0][0].legend(ncol=2,fontsize=img.xy_lb_size-6)
# img.axs[0][1].legend(ncol=1,fontsize=img.xy_lb_size-6)
# img.axs[0][2].legend(ncol=1,fontsize=img.xy_lb_size-6)
#
# img.axs[0][0].set_title("Fore Z:%.3f. Back $Z_1$: %.3f"%(z[fore_tag], z[len_tag1]),fontsize=img.xy_lb_size-6)
# img.axs[0][1].set_title("Zoom in",fontsize=img.xy_lb_size-6)
# img.axs[0][2].set_title("$10^2\Delta R$",fontsize=img.xy_lb_size-6)
#
#
# img.set_label(0, 0, 0, "$R = \\frac{1-\chi_L / \chi_S(Z_1)}{1-\chi_L / \chi_S(Z_2)}$")
# img.set_label(0, 0, 1, "$Z_2$")
# img.set_label(0, 1, 1, "$Z_2$")
# img.set_label(0, 2, 1, "$Z_2$")
#
# img.save_img(tool_box.file_name("E:/works/Galaxy-Galaxy_lensing/gamma_ratio/gamma_ratio.png"))
# img.show_img()
# img.close_img()
#
# img = Image_Plot(fig_x=8,fig_y=6)
# img.subplots(1, 1)
# for i in range(npara):
#     img.axs[0][0].plot(z,dist_c[i]/1000,c="C%d"%i,label="$\Omega_{m0}=%.2f$"%(Om0 + i * 0.02), linewidth=2.5)
# img.axs[0][0].legend(ncol=1,fontsize=img.xy_lb_size-6)
# img.axs[0][0].set_title("Comoving distance",fontsize=img.xy_lb_size-6)
# img.set_label(0, 0, 1, "$Z$")
# img.set_label(0, 0, 0, "$\chi \  \\rm{[10^3 Mpc]}$",size=img.xy_lb_size-4)
# img.save_img("E:/works/Galaxy-Galaxy_lensing/gamma_ratio/dist.png")
# img.show_img()
# img.close_img()
exit()


gamma = numpy.zeros((npara, num))
z_source = numpy.linspace(z_lens, z_max, num)

img = Image_Plot()
img.subplots(1, 1)

axins = zoomed_inset_axes(img.axs[0][0], 1.8, loc=7)

for i in range(npara):
    Om0_ = Om0 + i * 0.02
    cosmos = FlatLambdaCDM(H0=H0, Om0=Om0_)
    dist_lens = cosmos.comoving_distance(z_lens).value
    print(gamma_coeff(z_lens, delta_sigma, H0, Om0_))
    for j in range(num):
        dist = cosmos.comoving_distance(z_source[j]).value
        gamma[i, j] = gamma_coeff(z_lens, delta_sigma, H0, Om0_) * (1 - dist_lens / dist)*100
    img.axs[0][0].plot(z_source, gamma[i], label="Om0=%.2f" % Om0_)

    axins.plot(z_source, gamma[i])

axins.set_xlim(1.5,z_max)
axins.set_ylim(1.7,2.2)
mark_inset(img.axs[0][0], axins, loc1=2, loc2=4, fc="none", ec="black")
# ip = mark_inset.InsetPosition(img.axs[0][0], [0.5, 0.1, 0.4, 0.2])
# axins.set_axes_locator(ip)

img.axs[0][0].legend(ncol=2, fontsize=img.xy_lb_size - 7,bbox_to_anchor=(0.1, 0.25))
img.set_label(0, 0, 0, "$\gamma_t\\times10^2$")
img.set_label(0, 0, 1, "Z")
# img.axs[0][0].set_yscale("log")
img.axs[0][0].set_title("$z_{lens} = %.2f,\ H0 = 70,\  \Delta\Sigma=%.1f \ \\rm{h \cdot M_\odot \cdot pc^{-2}}$"
                        % (z_lens, delta_sigma), fontsize=img.xy_lb_size - 5)
img.save_img("E:/gamma.png")
img.show_img()













exit()
lensing_low = numpy.loadtxt("E:/Github/astrophy-research/galaxy-galaxy lensing/GGL_calculation/lensing_low/data.dat")

cata_path = "F:/works/gg_lensing/cmass/fourier/old_cata/Tomo/"
h5f = h5py.File(cata_path+ "w_1_1.hdf5","r")
redshift = h5f["/Z"].value
h5f.close()
print(redshift.min(), redshift.max())
fq = Fourier_Quad(12, 124)
filen = 9
result = numpy.zeros((6, filen))
for i in range(filen):
    h5f = h5py.File(cata_path + "radius_%d.hdf5"%i,"r")
    data = h5f["/pair_data_0"].value
    h5f.close()
    crit_integ = data[:,4]*554.682135528
    mgt = data[:,0]*crit_integ
    mgx = data[:,1]*crit_integ
    mnut = data[:,2]
    mnux = data[:,3]
    print(data[:,-1].min())
    print(mgt.shape[0])
    img = Image_Plot()
    img.subplots(1, 2)
    gt, gt_sig = fq.fmin_g_new(mgt, mnut, 8, left=-180, right=180, fig_ax=img.axs[0][0])[:2]
    gx, gx_sig = fq.fmin_g_new(mgx, mnux, 8, left=-180, right=180, fig_ax=img.axs[0][1])[:2]
    # print(result[:,i])
    img.show_img()
    result[0,i] = gt
    result[1,i] = gt_sig
    result[2,i] = gx
    result[3,i] = gx_sig
    result[4,i] = data[:,-2].mean()


img = Image_Plot()
img.subplots(1,1)

names = ["old cata", "new cata",  "lensing low"]

img.axs[0][0].errorbar(result[4], result[0], result[1],capsize=3, c="C1",label="T")
img.axs[0][0].errorbar(result[4], result[2], result[3],capsize=3, c="C2",label="x")
img.axs[0][0].errorbar(lensing_low[:,0], lensing_low[:,1], lensing_low[:,2], capsize=3,c="C4",label="%s"%names[2])
img.axs[0][0].legend()
img.axs[0][0].set_yscale("log")
img.axs[0][0].set_xscale("log")
img.axs[0][0].set_ylim(0.1,180)
img.show_img()