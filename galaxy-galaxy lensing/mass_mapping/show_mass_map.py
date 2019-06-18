import platform
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import matplotlib.pyplot as plt
from matplotlib import cm
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
path.append("E:/Github/astrophy-research/mylib/")
import plot_tool
import tool_box
import h5py


source_no = argv[1]
crit_z = float(argv[2])
file_name = argv[3]
data_path = "/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/result/mass_map/CFHT_cluster/w_1/filter_%s/source_%s/"%(file_name,source_no)
h5f = h5py.File("/mnt/ddnfs/data_users/hkli/CFHT/gg_lensing/data/redshift.hdf5","r")

# data_path = "G:\Backup\works\Mass_Mapping/source_2/".replace("\\","/")
# h5f = h5py.File("G:/Backup/works/Mass_Mapping/redshift.hdf5","r")


dist_refer = h5f["distance"].value
red_refer = h5f["redshift"].value
h5f.close()

gets_item = [["para", "nx", "0"], ["para", "ny", "0"],
             ["para", "delta RA (arcmin)", "0"],["para", "delta Dec (arcmin)", "0"]]
para_items = tool_box.config(data_path+"result.dat", ["get" for i in range(len(gets_item))], gets_item)

nx, ny, d_dx, d_dy = int(para_items[0]), int(para_items[1]), float(para_items[2]),float(para_items[3])

np_data = numpy.load(data_path + "result.npz")

block_num = 10

data_set = np_data["arr_0"]
sp = data_set.shape

redshift_bin_num = int(sp[0]/block_num/ny)

ra_bin = np_data["arr_1"]
dec_bin = np_data["arr_2"]
redshift_bin = np_data["arr_3"]
foregal = np_data["arr_4"]

inverse = range(ny-1,-1,-1)

sgima_coeff = 388.283351
tag = tool_box.find_near(red_refer, foregal[2])
dist_len = dist_refer[tag]
tag = tool_box.find_near(red_refer, crit_z)
dist_source = dist_refer[tag]

crit_density_z = sgima_coeff*dist_source/(dist_source - dist_len)/dist_len*(1+foregal[2])
print("Foregal Z: %f"%foregal[2])
print("\Sigma_c(z=%.2f) = %.4f"%(crit_z, crit_density_z))
print("%d redshift bins\n"%redshift_bin_num,redshift_bin)
if crit_z <= foregal[2]+0.1:
    print("Too small Z")
    exit()

for ir in range(redshift_bin_num):

    result = data_set[ir*block_num*ny:(ir+1)*block_num*ny]

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

    gamma_tc = dens_t/crit_density_z
    gamma_xc = dens_x/crit_density_z

    kappa_recon_ks95 = tool_box.shear2kappa_ks95(gamma_t,gamma_x).real
    kappa_recon_ks95_crit = tool_box.shear2kappa_ks95(gamma_tc, gamma_xc).real

    kappa_recon = tool_box.shear2kappa(gamma_t, gamma_x)
    kappa_recon_crit = tool_box.shear2kappa(gamma_tc, gamma_xc)

    idx_t1 = numpy.abs(gamma_t) < 1.e-9
    idx_t2 = gamma_t < -1
    idx_t3 = gamma_t_sig < -1
    idx_x1 = numpy.abs(gamma_x) < 1.e-9
    idx_x2 = gamma_x < -1
    idx_x3 = gamma_x_sig < -1
    idx_n = num < 1

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
    gamma_c_sig = numpy.sqrt(dens_t_sig**2/crit_density_z**2 + dens_x_sig**2*dens_x**2/crit_density_z**4)

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

    density_data = [[dens_t, dens_x, sigma_mean], [dens_t_sig, dens_x_sig, num]]
    dens_titles = [["$< g_1\Sigma_c>$", "$< g_2\Sigma_c>$", "$<\Sigma_c>$"],
             ["$\delta <g_1\Sigma_c>$","$\delta <g_2\Sigma_c>$","$number$"]]

    gamma_data = [[gamma_t, gamma_x, gamma], [gamma_t_sig, gamma_x_sig, num],[gamma_tc, gamma_xc, num]]
    gamma_titles = [["$g_1$", "$g_2$", "$g$"],
              ["$\delta g_1$", "$\delta g_2$", "$number$"],
              ["$< g_1\Sigma_c>/\Sigma_{z=%.2f}$"%crit_z,"$< g_2\Sigma_c>/\Sigma_{z=%.2f}$"%crit_z,"nothing"]]

    kappa_data = [[kappa_recon_ks95, kappa_recon_ks95_crit], [kappa_recon, kappa_recon_crit]]
    kappa_titles = [["$\kappa KS\_95 Fourier$ ", "$\kappa KS\_95 Fourier critical$"],
              ["$\kappa real$", "$\kappa real criticak$"]]

    # show density map
    img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
    img.subplots(2,3)
    cmap = plt.get_cmap('YlOrRd')
    sm = plt.cm.ScalarMappable(cmap=cmap)

    for i in range(2):
        for j in range(3):
            density_data[i][j][idx] = 0
            ax = img.axs[i][j].imshow(density_data[i][j][inverse],cmap=cmap)
            img.set_label(i,j, 1, "RA")
            img.set_label(i,j, 0, "DEC")
            img.axs[i][j].set_title(dens_titles[i][j],fontsize=img.xy_lb_size)
            img.figure.colorbar(ax, ax=img.axs[i][j])

    img.save_img(data_path + "density_%d.png"%ir)
    if platform.system() != 'Linux':
        img.show_img()
    img.close_img()

    # show shear map
    img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
    img.subplots(3, 3)
    cmap = plt.get_cmap('YlOrRd')
    sm = plt.cm.ScalarMappable(cmap=cmap)

    for i in range(3):
        for j in range(3):
            gamma_data[i][j][idx] = 0
            print(mask[idx])
            ax = img.axs[i][j].imshow(gamma_data[i][j][inverse], cmap=cmap)
            img.set_label(i, j, 1, "RA")
            img.set_label(i, j, 0, "DEC")
            img.axs[i][j].set_title(gamma_titles[i][j], fontsize=img.xy_lb_size)
            img.figure.colorbar(ax, ax=img.axs[i][j])

    img.save_img(data_path + "shear_%d.png"%ir)
    if platform.system() != 'Linux':
        img.show_img()
    img.close_img()

    # show the recovered kappa map
    img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
    img.subplots(2,2)

    for i in range(2):
        for j in range(2):
            ax = img.axs[i][j].imshow(kappa_data[i][j][inverse],cmap=cm.jet)
            img.set_label(i, j, 1, "RA")
            img.set_label(i, j, 0, "DEC")
            img.axs[i][j].set_title(kappa_titles[i][j], fontsize=img.xy_lb_size)
            img.figure.colorbar(ax, ax=img.axs[i][j])

    img.save_img(data_path + "kappa_%d.png"%ir)
    if platform.system() != 'Linux':
        img.show_img()
    img.close_img()

    max_g = gamma[numpy.abs(gamma)<0.1].max()
    max_gc = gamma_c[numpy.abs(gamma_c) < 0.1].max()
    print(max_g,max_gc)
    dec_sq_len = (dec_bin[2] - dec_bin[1])
    ra_sq_len = (ra_bin[2] - ra_bin[1])
    max_len = ra_sq_len*0.9

    dg_scale = gamma/max_g*max_len/2
    dg_scale_c = gamma_c / max_gc * max_len / 2

    # plot the shear field
    img = plot_tool.Image_Plot(fig_x=20, fig_y=20)
    img.subplots(1,2)

    shear_bench = 0.03
    scale_len = shear_bench/max_g*max_len

    x1, x2 = ra_min + ra_sq_len*6, ra_min + ra_sq_len*6 + scale_len
    y1,y2 = dec_max + dec_sq_len*3, dec_max + dec_sq_len*3
    img.axs[0][0].plot([x1,x2],[y1,y2],c="black")
    img.axs[0][0].text(0.03, 0.93, "shear=%.3f"%shear_bench, color='black', ha='left',
            va='center', transform=img.axs[0][0].transAxes, fontsize=img.legend_size-5)

    img.axs[0][0].scatter(foregal[0], foregal[1],s=200,facecolors="none",edgecolors="r",marker="*")
    img.axs[0][1].scatter(foregal[0], foregal[1], s=200, facecolors="none", edgecolors="r", marker="*")
    for i in range(ny + 1):
        img.axs[0][0].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--",alpha=0.5,linewidth=0.3)
        img.axs[0][1].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--",
                           alpha=0.5, linewidth=0.3)
    for j in range(nx + 1):
        img.axs[0][0].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black",linestyle="--" ,alpha=0.5, linewidth=0.3)
        img.axs[0][1].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black", linestyle="--",
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

                x = (ra_bin[j] + ra_bin[j+1])/2
                y = (dec_bin[i] + dec_bin[i+1])/2
                if color == 1:
                    img.axs[0][0].plot([x+dx, x-dx], [y+dy, y-dy],c=cl)
                    img.axs[0][1].plot([x + dxc, x - dxc], [y + dyc, y - dyc], c=cl)
                else:
                    img.axs[0][0].plot([x + dx, x - dx], [y + dy, y - dy], c="C0")
                    img.axs[0][1].plot([x + dxc, x - dxc], [y + dyc, y - dyc], c="C0")
    img.axs[0][0].set_title("From the estimated $g$",fontsize=img.xy_lb_size)
    img.axs[0][1].set_title("Recoved from the $< g_{1/2}\Sigma_c>/\Sigma_{z=%.2f}$"%crit_z,fontsize=img.xy_lb_size)
    img.set_label(0, 0, 1, "RA")
    img.set_label(0, 0, 0, "DEC")
    img.set_label(0, 1, 1, "RA")
    img.set_label(0, 1, 0, "DEC")
    if color == 1:
        sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
        sm._A = []
        plt.colorbar(sm, ax=img.axs[0][0])

    img.save_img(data_path + "shear_map_%d_z=%.2f.png"%(ir,crit_z))
    if platform.system() != 'Linux':
        img.show_img()