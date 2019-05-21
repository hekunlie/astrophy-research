import numpy
import matplotlib.pyplot as plt
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
import plot_tool
import tool_box



data_path = "E:/mass_map/CFHT_cluster/w_1/source_0/"

gets_item = [["para", "nx", "0"], ["para", "ny", "0"]]
para_items = tool_box.config(data_path+"result.dat", ["get", "get"], gets_item)

nx, ny = int(para_items[0]), int(para_items[1])

np_data = numpy.load(data_path + "result.npz")

result = np_data["arr_0"]
ra_bin = np_data["arr_1"]
dec_bin = np_data["arr_2"]
foregal = np_data["arr_3"]

inverse = range(ny-1,-1,-1)

gamma_t = result[:ny]
gamma_t_sig = result[ny:2*ny]
gamma_x = result[2*ny:3*ny]
gamma_x_sig = result[3*ny:4*ny]
num = result[4*ny:5*ny]
position_angle = result[5*ny:6*ny]


# g1 = gamma_t*numpy.cos(-2*position_angle) - gamma_x*numpy.sin(-2*position_angle)
# g2 = gamma_t*numpy.sin(-2*position_angle) - gamma_x*numpy.cos(-2*position_angle)
gamma = numpy.sqrt(gamma_t**2 + gamma_x**2)
gamma_sig = numpy.abs(gamma_t)/gamma*gamma_t_sig + numpy.abs(gamma_x)/gamma*gamma_x_sig

img_theta = plot_tool.Image_Plot(fig_x=8, fig_y=8)
img_theta.create_subfig(1, 1)
# position angle
theta = numpy.arccos(numpy.sqrt((1-gamma_t/gamma)/2))

ax = img_theta.axs[0][0].imshow(theta)
img_theta.figure.colorbar(ax, ax=img_theta.axs[0][0])
img_theta.show_img()
plt.close()

ra_min, ra_max = ra_bin.min(), ra_bin.max()
dec_min, dec_max = dec_bin.min(), dec_bin.max()

datas = [[gamma_t, gamma_x, position_angle/numpy.pi*180], [gamma_t_sig, gamma_x_sig, num]]
titles = [["$\gamma_t$", "$\gamma_x$", "$\\theta$"],
         ["$\delta \gamma_t$","$\delta \gamma_x$","$number$"]]
img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
img.create_subfig(2,3)
cmap = plt.get_cmap('YlOrRd')
sm = plt.cm.ScalarMappable(cmap=cmap)

for i in range(2):
    for j in range(3):
        ax = img.axs[i][j].imshow(datas[i][j][inverse],cmap=cmap)
        img.tick_label(i,j, 1, "RA")
        img.tick_label(i,j, 0, "DEC")
        img.axs[i][j].set_title(titles[i][j],fontsize=img.xy_lb_size)
        img.figure.colorbar(ax, ax=img.axs[i][j])

img.save_img(data_path + "result.png")
img.show_img()
plt.close()

img = plot_tool.Image_Plot(fig_x=12, fig_y=10)
img.create_subfig(1,2)

img.axs[0][0].scatter(foregal[0], foregal[1],s=200,facecolors="none",edgecolors="r",marker="*")
for i in range(ny + 1):
    img.axs[0][0].plot([ra_min, ra_max], [dec_bin[i], dec_bin[i]], c="black", linestyle="--",alpha=0.5,linewidth=0.3)
for j in range(nx + 1):
    img.axs[0][0].plot([ra_bin[j], ra_bin[j]], [dec_min, dec_max], c="black",linestyle="--" ,alpha=0.5, linewidth=0.3)

# gamma[4,4] = 0

ax = img.axs[0][1].imshow(gamma[inverse], cmap=cmap)
img.tick_label(0, 1, 1, "RA")
img.tick_label(0, 1, 0, "DEC")
img.figure.colorbar(ax, ax=img.axs[0][1])

max_g = gamma.max()
max_len = (ra_bin[2] - ra_bin[1])*0.7

dg_scale = gamma/max_g*max_len/2

new_sig = numpy.zeros_like(gamma_sig)
new_sig[0:ny, 0:nx] = gamma_sig
# new_sig[4,4] = 0.01

norm = plt.Normalize(vmin=numpy.min(new_sig), vmax=numpy.max(new_sig))
cmap = plt.get_cmap('plasma')

color = 2
for i in range(ny):
    for j in range(nx):
        # if i == 4 and j == 4:
        #     pass
        # else:
        dx = numpy.abs(dg_scale[i,j]*numpy.cos(theta[i,j]))
        dy = numpy.abs(dg_scale[i,j]*numpy.sin(theta[i,j]))
        cl = cmap(norm(new_sig[i,j]))

        x = (ra_bin[j] + ra_bin[j+1])/2
        y = (dec_bin[i] + dec_bin[i+1])/2
        if color == 1:
            img.axs[0][0].plot([x+dx, x-dx], [y+dy, y-dy],c=cl)
        else:
            img.axs[0][0].plot([x + dx, x - dx], [y + dy, y - dy], c="C0")
        img.tick_label(0, 0, 1, "RA")
        img.tick_label(0, 0, 0, "DEC")
if color == 1:
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm._A = []
    plt.colorbar(sm, ax=img.axs[0][0])

img.save_img(data_path + "shear_map.png")
img.show_img()