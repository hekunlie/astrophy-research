import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import matplotlib.pyplot as plt
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
import plot_tool
import tool_box
from scipy import signal


data_path = "G:\Backup\works\Mass_Mapping/source_2/".replace("\\","/")

gets_item = [["para", "nx", "0"], ["para", "ny", "0"],
             ["para", "delta RA (arcmin)", "0"],["para", "delta Dec (arcmin)", "0"]]
para_items = tool_box.config(data_path+"result.dat", ["get" for i in range(len(gets_item))], gets_item)

nx, ny, d_dx, d_dy = int(para_items[0]), int(para_items[1]), float(para_items[2]),float(para_items[3])

np_data = numpy.load(data_path + "result.npz")

block_num = 10

result = np_data["arr_0"]
sp = result.shape

redshift_bin_num = int(sp[0]/block_num/ny)

ra_bin = np_data["arr_1"]
dec_bin = np_data["arr_2"]
redshift_bin = np_data["arr_3"]
foregal = np_data["arr_4"]

print("Foregal Z: %f"%foregal[2])
print("%d redshift bins\n"%redshift_bin_num,redshift_bin)

inverse = range(ny-1,-1,-1)

for ir in range(redshift_bin_num):

    sub_data = result[ir*block_num*ny:(ir+1)*block_num*ny]

    dens_t = result[:ny]
    dens_t_sig = result[ny:2*ny]
    dens_x = result[2*ny:3*ny]
    dens_x_sig = result[3*ny:4*ny]
    sigma_mean = result[4*ny:5*ny]

    gamma_t = result[5*ny:6*ny]
    gamma_t_sig = result[6*ny:7*ny]
    gamma_x = result[7*ny:8*ny]
    gamma_x_sig = result[8*ny:9*ny]
    num = result[9*ny:10*ny]

    gamma_tc = dens_t/sigma_mean
    gamma_xc = dens_x/sigma_mean

    idx_t1 = numpy.abs(gamma_t) < 1.e-9
    idx_t2 = gamma_t < -1
    idx_t3 = gamma_t_sig < -1
    idx_x1 = numpy.abs(gamma_x) < 1.e-9
    idx_x2 = gamma_x < -1
    idx_x3 = gamma_x_sig < -1
    idx_n = num < 0

    idxs = [idx_t1,idx_t2,idx_t3,idx_x1,idx_x2,idx_x3, idx_n]

    mask = numpy.arange(0, ny*nx).reshape((ny,nx))

    idx = idx_t1 + idx_t2 + idx_t3 + idx_x1 + idx_x2 + idx_x3 + idx_n
    idx_c = idx_t1 & idx_t2 & idx_t3 & idx_x1 & idx_x2 & idx_x3
    if idx_c.sum() > 0 or idx.sum() > 0:
        print(idx.sum(),idx_c.sum())
        continue

    gamma = numpy.sqrt(gamma_t**2 + gamma_x**2)
    gamma_sig = numpy.sqrt((gamma_t/gamma)**2*gamma_t_sig**2 + (gamma_x/gamma)**2*gamma_x_sig**2)
    gamma_c = numpy.sqrt(gamma_tc**2 + gamma_xc**2)
    gamma_c_sig = numpy.sqrt(dens_t_sig**2/sigma_mean**2 + dens_x_sig**2*dens_x**2/sigma_mean**4)

    # position angle
    cos_theta = numpy.sqrt((1+gamma_t/gamma)/2)
    sin_theta = gamma_x/2/gamma/cos_theta

    cos_theta_c = numpy.sqrt((1+gamma_tc/gamma_c)/2)
    sin_theta_c = gamma_xc/2/gamma_c/cos_theta_c

    theta = numpy.arccos(cos_theta)

    # gamma_t[idx] = 0
    # gamma_t_sig[idx] = 0
    # gamma_x[idx] = 0
    # gamma_x_sig[idx] = 0
    # # theta[idx] = 0
    # num[idx] = 0

    ra_min, ra_max = ra_bin.min(), ra_bin.max()
    dec_min, dec_max = dec_bin.min(), dec_bin.max()

    datas = [[dens_t, dens_x, sigma_mean], [dens_t_sig, dens_x_sig, num]]
    titles = [["$< g_1\Sigma_c>$", "$< g_2\Sigma_c>$", "$<\Sigma_c>$"],
             ["$\delta <g_1\Sigma_c>$","$\delta <g_2\Sigma_c>$","$number$"]]
    img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
    img.create_subfig(2,3)
    cmap = plt.get_cmap('YlOrRd')
    sm = plt.cm.ScalarMappable(cmap=cmap)

    for i in range(2):
        for j in range(3):
            datas[i][j][idx] = 0
            ax = img.axs[i][j].imshow(datas[i][j][inverse],cmap=cmap)
            img.tick_label(i,j, 1, "RA")
            img.tick_label(i,j, 0, "DEC")
            img.axs[i][j].set_title(titles[i][j],fontsize=img.xy_lb_size)
            img.figure.colorbar(ax, ax=img.axs[i][j])

    img.save_img(data_path + "density_%d.png"%ir)
    img.show_img()
    plt.close()

    datas = [[gamma_t, gamma_x, gamma], [gamma_t_sig, gamma_x_sig, num],[gamma_tc, gamma_xc, num]]
    titles = [["$g_1$", "$g_2$", "$g$"],
              ["$\delta g_1$", "$\delta g_2$", "$number$"],
              ["$< g_1\Sigma_c>/<\Sigma_c>$","$< g_2\Sigma_c>/<\Sigma_c>$","nothing"]]
    img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
    img.create_subfig(3, 3)
    cmap = plt.get_cmap('YlOrRd')
    sm = plt.cm.ScalarMappable(cmap=cmap)

    for i in range(3):
        for j in range(3):
            datas[i][j][idx] = 0
            print(mask[idx])
            ax = img.axs[i][j].imshow(datas[i][j][inverse], cmap=cmap)
            img.tick_label(i, j, 1, "RA")
            img.tick_label(i, j, 0, "DEC")
            img.axs[i][j].set_title(titles[i][j], fontsize=img.xy_lb_size)
            img.figure.colorbar(ax, ax=img.axs[i][j])

    img.save_img(data_path + "shear_%d.png"%ir)
    img.show_img()
    plt.close()


    max_g = gamma[numpy.abs(gamma)<0.1].max()
    max_gc = gamma_c[numpy.abs(gamma_c) < 0.1].max()
    print(max_g,max_gc)
    max_len = (ra_bin[2] - ra_bin[1])*60*0.7

    dg_scale = gamma/max_g*max_len/2
    dg_scale_c = gamma_c / max_gc * max_len / 2

    img = plot_tool.Image_Plot(fig_x=14, fig_y=14)
    img.create_subfig(2,1)

    img.axs[0][0].scatter(foregal[0]*60, foregal[1]*60,s=200,facecolors="none",edgecolors="r",marker="*")
    img.axs[1][0].scatter(foregal[0] * 60, foregal[1] * 60, s=200, facecolors="none", edgecolors="r", marker="*")
    for i in range(ny + 1):
        img.axs[0][0].plot([ra_min*60, ra_max*60], [dec_bin[i]*60, dec_bin[i]*60], c="black", linestyle="--",alpha=0.5,linewidth=0.3)
        img.axs[1][0].plot([ra_min * 60, ra_max * 60], [dec_bin[i] * 60, dec_bin[i] * 60], c="black", linestyle="--",
                           alpha=0.5, linewidth=0.3)
    for j in range(nx + 1):
        img.axs[0][0].plot([ra_bin[j]*60, ra_bin[j]*60], [dec_min*60, dec_max*60], c="black",linestyle="--" ,alpha=0.5, linewidth=0.3)
        img.axs[1][0].plot([ra_bin[j] * 60, ra_bin[j] * 60], [dec_min * 60, dec_max * 60], c="black", linestyle="--",
                           alpha=0.5, linewidth=0.3)
    norm = plt.Normalize(vmin=numpy.min(max_g), vmax=numpy.max(max_g))
    cmap = plt.get_cmap('plasma')

    print("abnormal: ",idx.sum())
    color = 2
    for i in range(ny):
        for j in range(nx):
            if not idx[i,j]:

                if gamma_x[i,j] < 0:
                    dx = -dg_scale[i,j]*cos_theta[i,j]
                    dy = -dg_scale[i,j]*sin_theta[i,j]
                else:
                    dx = dg_scale[i, j] * cos_theta[i, j]
                    dy = dg_scale[i, j] * sin_theta[i, j]

                if gamma_xc[i,j] < 0:
                    dxc = -dg_scale_c[i,j]*cos_theta_c[i,j]
                    dyc = -dg_scale_c[i,j]*sin_theta_c[i,j]
                else:
                    dxc = dg_scale_c[i, j] * cos_theta_c[i, j]
                    dyc = dg_scale_c[i, j] * sin_theta_c[i, j]
                cl = cmap(norm(gamma[i,j]))

                x = (ra_bin[j] + ra_bin[j+1])/2*60
                y = (dec_bin[i] + dec_bin[i+1])/2*60
                if color == 1:
                    img.axs[0][0].plot([x+dx, x-dx], [y+dy, y-dy],c=cl)
                    img.axs[1][0].plot([x + dxc, x - dxc], [y + dyc, y - dyc], c=cl)
                else:
                    img.axs[0][0].plot([x + dxc, x - dxc], [y + dyc, y - dyc], c="C0")
                    img.axs[1][0].plot([x + dxc, x - dxc], [y + dyc, y - dyc], c="C0")

    img.tick_label(0, 0, 1, "RA")
    img.tick_label(1, 0, 0, "DEC")
    img.tick_label(1, 0, 1, "RA")
    img.tick_label(1, 0, 0, "DEC")
    if color == 1:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        plt.colorbar(sm, ax=img.axs[0][0])

    img.save_img(data_path + "shear_map_%d.png"%ir)
    img.show_img()