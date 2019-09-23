from sys import path
path.append("/home/hklee/work/mylib/")
path.append("D:/GitHub/astrophy-research/mylib/")
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import h5py
import numpy
import gauss_fit_fun


coeff = 554.682135528

fq = Fourier_Quad(12,123)

total_path = "E:/GGL_FQ/"
# signal = numpy.zeros((2,10))

npz_cache = numpy.load("%scache.npz"%total_path)
signal = npz_cache["arr_0"]
for ig in range(10):
    h5f = h5py.File("%sradius_%d.hdf5"%(total_path,ig),"r")
    print(list(h5f.keys()))
    data_col = h5f["/data_col"].value
    pair_count = h5f["/pair_count"].value
    data = h5f["/pair_data_0"].value
    radius_bin = h5f["/radius_bin"].value
    h5f.close()


    # print(data_col)
    # print(pair_count)
    # print(data.shape)
    # print(radius_bin.shape)
    # print(data[:5])

    sigma = data[:,4]*coeff
    mg1 = data[:,0]
    mg2 = data[:,1]
    mnu1 = data[:,2]
    mnu2 = data[:,3]
    radius_sep = data[:, 5]
    redshift = data[:, 6]

    # gh, gh_sig = fq.find_shear(mg1*sigma, mnu1, 8, left=0, right=120)[:2]
    # print(gh, gh_sig)

    img = Image_Plot()
    img.subplots(2, 3)
    img.set_style()
    img.axs[0][0].scatter(mg1, radius_sep, s=3)
    img.axs[0][1].scatter(mg1, redshift, s=3)
    img.axs[0][2].scatter(mg1, sigma, s=3)
    #
    # bin_num = 51
    # bin_x = numpy.linspace(mg1.min(),mg1.max(), bin_num)
    # bin_y = numpy.linspace(radius_sep.min(),radius_sep.max(),bin_num)
    # num_bin = numpy.histogram2d(radius_sep, mg1, bins=[bin_y, bin_x])[0]

    fig = img.axs[1][0].hist2d(mg1, radius_sep, bins=[30,30])[3]
    img.figure.colorbar(fig, ax=img.axs[1][0])

    fig = img.axs[1][1].hist2d(mg1, redshift, bins=[30,30])[3]
    img.figure.colorbar(fig, ax=img.axs[1][1])

    fig = img.axs[1][2].hist2d(mg1,sigma,bins=[30,30])[3]
    img.figure.colorbar(fig, ax=img.axs[1][2])

    for i in range(2):
        img.set_label(i, 0, 0, "Separation", size=img.xy_lb_size)
        img.set_label(i, 1, 0, "Z", size=img.xy_lb_size)
        img.set_label(i, 2, 0, "$D_S / (D_LD_{LS})$", size=img.xy_lb_size)
        for j in range(3):
            img.set_label(i,j,1,"$\hat G$", size=img.xy_lb_size)

    img.save_img(total_path + "2d_hist_%d.png"%ig)
    # img.show_img()

    continue
    bins = numpy.histogram(mg1, 3000)[1]
    g_guess = numpy.linspace(-0.001, 0.2, 15 + i*5) * coeff
    num_move = gauss_fit_fun.get_flow(mg1, mnu1, g_guess)

    img = Image_Plot()
    img.subplots(1,1)
    img.set_style()
    # for j in range(len(g_guess)):
    #     img.axs[0][0].hist(mg1 - g_guess[j]*mnu1, bins, histtype="step", label="%f"%g_guess[j])
    # # img.axs[0][0].set_xscale("symlog")
    # img.axs[0][0].set_xlim(-500,500)
    # img.axs[0][0].legend()

    img.axs[0][0].scatter(g_guess[1:], num_move)
    img.axs[0][0].set_title("$\Delta \Sigma = %.3f (%.3f)$"%(gh, gh_sig))
    img.set_label(0,0,0,"NUM_CROSS_ZERO")
    img.set_label(0,0,1,"$\hat{\Delta \Sigma}$")
    img.save_img(total_path + "move_%d.png"%i)
    # img.show_img()


    # signal[0,i] = gh*coeff
    # signal[1,i] = gh_sig*coeff

    # print(signal[:,i])

# numpy.savez("%scache.npz"%total_path, signal, radius_bin)
#
# img = Image_Plot()
# img.subplots(1,1)
# img.axs[0][0].errorbar(radius_bin[:10,0], signal[0], signal[1])
# img.axs[0][0].set_xscale("log")
# img.axs[0][0].set_yscale("log")
# img.show_img()