import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
import matplotlib.pyplot as plt
from sys import path
path.append("E:/Github/astrophy-research/mylib/")
import plot_tool
import tool_box
from scipy import signal

def get_kappa(g1, g2, dx, dy):

    ny, nx = g1.shape

    g = g1 ** 2 + g2 ** 2
    g1_dx = (g1[0:ny - 1, 1:nx] - g1[0:ny - 1, 0:nx - 1]) / dx
    g1_dy = (g1[1:ny, :nx - 1] - g1[0:ny - 1, :nx - 1]) / dy

    g2_dx = (g2[0:ny - 1, 1:nx] - g2[0:ny - 1, 0:nx - 1]) / dx
    g2_dy = (g2[1:ny, :nx - 1] - g2[0:ny - 1, :nx - 1]) / dy

    coeff_0 = 1 - g[:ny - 1, :nx - 1]
    coeff_11 = (1 + g1[:ny - 1, :nx - 1]) / coeff_0
    coeff_12 = g2[:ny - 1, :nx - 1] / coeff_0
    coeff_21 = coeff_12
    coeff_22 = (1 - g1[:ny - 1, :nx - 1]) / coeff_0

    df_dx = (coeff_11 * (g1_dx + g2_dy) + coeff_12 * (g2_dx - g1_dy))*dx
    df_dy = (coeff_21 * (g1_dx + g2_dy) + coeff_22 * (g2_dx - g1_dy))*dy
    print(df_dx.shape, df_dy.shape)
    margin = 2
    kappa = numpy.zeros((ny + margin, nx + margin))

    fx = numpy.zeros((ny-1,nx-1))
    fy = numpy.zeros_like(fx)
    for i in range(nx-1):
        fx[:,i] = df_dx[:,:i+1].sum(axis=1)
    for i in range(ny-1):
        fy[i] = df_dy[:i+1].sum(axis=0)
    kappa[1:ny,1:nx] = 1 - numpy.exp(fx+fy)

    ld = 10
    my, mx = numpy.mgrid[0:ld, 0:ld]
    w = 5
    gauss = 1./2/numpy.pi/w/w*numpy.exp(-((my-5.5)**2+(mx-5.5)**2)/2/w**2)
    kappa_s = signal.convolve(kappa,gauss,mode="same")

    ky, kx = kappa.shape
    inverse = range(ky - 1, -1, -1)
    plt.subplot(121)
    plt.imshow(kappa[inverse])
    plt.subplot(122)
    plt.imshow(kappa_s[inverse])
    plt.show()
    plt.close()

    return df_dx, df_dy

data_path = "F:/works/w_1/source_2/"

gets_item = [["para", "nx", "0"], ["para", "ny", "0"],
             ["para", "delta RA (arcmin)", "0"],["para", "delta Dec (arcmin)", "0"]]
para_items = tool_box.config(data_path+"result.dat", ["get" for i in range(len(gets_item))], gets_item)

nx, ny, d_dx, d_dy = int(para_items[0]), int(para_items[1]), float(para_items[2]),float(para_items[3])

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


gamma = numpy.sqrt(gamma_t**2 + gamma_x**2)
gamma_sig = numpy.abs(gamma_t)/gamma*gamma_t_sig + numpy.abs(gamma_x)/gamma*gamma_x_sig

gamma_sig[idx] = 0
#
#
# position angle
cos_theta = numpy.sqrt((1+gamma_t/gamma)/2)
sin_theta = gamma_x/2/gamma/cos_theta

theta = numpy.arccos(cos_theta)

# gamma_t[idx] = 0
# gamma_t_sig[idx] = 0
# gamma_x[idx] = 0
# gamma_x_sig[idx] = 0
# # theta[idx] = 0
# num[idx] = 0

ra_min, ra_max = ra_bin.min(), ra_bin.max()
dec_min, dec_max = dec_bin.min(), dec_bin.max()

datas = [[gamma_t, gamma_x, num], [gamma_t_sig, gamma_x_sig, num]]
titles = [["$\gamma_t$", "$\gamma_x$", "$position angle$"],
         ["$\delta \gamma_t$","$\delta \gamma_x$","$number$"]]
img = plot_tool.Image_Plot(fig_x=8, fig_y=8)
img.create_subfig(2,3)
cmap = plt.get_cmap('YlOrRd')
sm = plt.cm.ScalarMappable(cmap=cmap)

for i in range(2):
    for j in range(3):
        datas[i][j][idx] = 0
        print(mask[idx])
        ax = img.axs[i][j].imshow(datas[i][j][inverse],cmap=cmap)
        img.tick_label(i,j, 1, "RA")
        img.tick_label(i,j, 0, "DEC")
        img.axs[i][j].set_title(titles[i][j],fontsize=img.xy_lb_size)
        img.figure.colorbar(ax, ax=img.axs[i][j])

img.save_img(data_path + "result.png")
img.show_img()
plt.close()

img = plot_tool.Image_Plot(fig_x=12, fig_y=12)
img.create_subfig(1,1)
gamma[idx] = 0
ax = img.axs[0][0].imshow(gamma[inverse], cmap=cmap)
img.tick_label(0, 0, 1, "RA")
img.tick_label(0, 0, 0, "DEC")
img.figure.colorbar(ax, ax=img.axs[0][0])
img.save_img(data_path + "gamma.png")
img.show_img()

max_g = gamma[numpy.abs(gamma)<0.1].max()
print(max_g)
max_len = (ra_bin[2] - ra_bin[1])*60*0.7

dg_scale = gamma/max_g*max_len/2

new_sig = numpy.zeros_like(gamma_sig)
new_sig[0:ny, 0:nx] = gamma_sig
# new_sig[4,4] = 0.01



get_kappa(gamma_t, gamma_x, d_dx, d_dy)


img = plot_tool.Image_Plot(fig_x=16, fig_y=16)
img.create_subfig(1,1)

img.axs[0][0].scatter(foregal[0]*60, foregal[1]*60,s=200,facecolors="none",edgecolors="r",marker="*")
for i in range(ny + 1):
    img.axs[0][0].plot([ra_min*60, ra_max*60], [dec_bin[i]*60, dec_bin[i]*60], c="black", linestyle="--",alpha=0.5,linewidth=0.3)
for j in range(nx + 1):
    img.axs[0][0].plot([ra_bin[j]*60, ra_bin[j]*60], [dec_min*60, dec_max*60], c="black",linestyle="--" ,alpha=0.5, linewidth=0.3)
norm = plt.Normalize(vmin=numpy.min(new_sig), vmax=numpy.max(new_sig))
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
            cl = cmap(norm(new_sig[i,j]))

            x = (ra_bin[j] + ra_bin[j+1])/2*60
            y = (dec_bin[i] + dec_bin[i+1])/2*60
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

img.save_img(data_path + "shear_map.pdf")
img.show_img()