import matplotlib
matplotlib.use("Agg")
import numpy
import os
my_home = os.popen("echo $HOME").readlines()[0][:-1]
from sys import path,argv
path.append('%s/work/mylib/'%my_home)
from Fourier_Quad import Fourier_Quad
from plot_tool import Image_Plot
import emcee
import corner
import time
import matplotlib.pyplot as plt
from multiprocessing import Pool
import h5py
from scipy.optimize import least_squares
import MCMC_program
import tool_box


expo_num = int(argv[1])
ncpus = int(argv[2])
shear_cmd = argv[3]


parent_path = "/mnt/perc/hklee/CFHT/multi_shear/cluster_field/"
# parent_path = "/mnt/ddnfs/data_users/hkli/CFHT/multi_shear/cluster_field/result/"

envs_path = parent_path + "param_slope.dat"
contents = [['param', "RA", "1"], ['param',"DEC", "1"],
            ['param', "a1", "1"], ['param', "a2", "1"], ['param', "a3", "1"]]

var_items = tool_box.config(envs_path, ['get' for i in range(len(contents))], contents)

# the field size in unit of arcmin
delta_ra = float(var_items[0])
delta_dec = float(var_items[1])
parameters = [float(var_items[2]), float(var_items[3]), float(var_items[4])]
print("PARAM:",parameters)
nwalkers = 300
ndim = len(parameters)
step = 1000
print("Walker: %d. Step: %d."%(nwalkers, step))
fq = Fourier_Quad(10, 112)
bin_num = 8
bin_num2 = int(bin_num / 2)

nx, ny = 100, 100
half_side = nx/2
pixel_scale = delta_ra/nx



for i in range(expo_num):
    h5f = h5py.File(parent_path + "result/expo_%d_slope.hdf5"%i,"r")
    temp = h5f["/data"].value
    h5f.close()
    if i == 0:
        data = temp
    else:
        data = numpy.row_stack((data, temp))


# radius = numpy.sqrt(data[:, 8]**2 + data[:, 9]**2)
# idx = radius <= fit_radius
# data = data[idx]

mg1 = data[:, 0]
mg2 = data[:, 1]
mnu1 = data[:, 2] + data[:, 3]
mnu2 = data[:, 2] - data[:, 3]
g1 = data[:, 5]
shear = g1
ra = data[:, 6]
dec = data[:, 7]

radius = numpy.sqrt(ra**2 + dec**2)

mg_bins = fq.set_bin(mg1, bin_num, 5)
inverse = range(bin_num2 - 1, -1, -1)
print("Bins:", mg_bins)

p0 = numpy.zeros((nwalkers, ndim))
p0[:, 0] = numpy.random.uniform(-0.01, 0.01, nwalkers)
p0[:, 1] = numpy.random.uniform(-0.01, 0.01, nwalkers)
p0[:, 2] = numpy.random.uniform(-0.01, 0.01, nwalkers)

if shear_cmd == "g1":
    shear_tag = 1
    MG = mg1
    NU = mnu1
elif shear_cmd == "g2":
    shear_tag = 2
    MG = mg2
    NU = mnu2
else:
    exit()
t1 = time.time()
with Pool(ncpus) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, MCMC_program.ln_prob_g_slope,
                                    args=[MG, NU, mg_bins, bin_num2, inverse, ra, dec],pool=pool)
    pos, prob, state = sampler.run_mcmc(p0, step)
t2 = time.time()
print("Time:  %2.f sec"%(t2 - t1))

# Plot the walkers
img = Image_Plot(fig_x=16, fig_y=4)
img.subplots(ndim,1)
for i in range(nwalkers):
    for j in range(ndim):
        img.axs[j][0].plot(range(step),sampler.chain[i, :, j], color='grey',alpha=0.6)
        img.axs[j][0].plot([0,step], [parameters[j], parameters[j]])
img.save_img(parent_path + "pic/%s_mcmc_walkers_nw_%d_stp_%d_expo_%d_slope.png"%(shear_cmd, nwalkers, step, expo_num))
img.close_img()

# Plot the panels
samples = sampler.chain[:, 150:, :].reshape((-1, ndim))
print(samples.shape)
corner_fig = plt.figure(figsize=(10, 10))
fig = corner.corner(samples, labels=["$a1$", "$a2$", "$a3$"], truths=parameters,
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".4f", title_kwargs={"fontsize": 12})
fig.savefig(parent_path + "pic/%s_mcmc_panel_nw_%d_stp_%d_expo_%d_slope.png"%(shear_cmd,nwalkers, step,expo_num))

fit_params = []
pr = numpy.percentile(samples, [16, 50, 84], axis=0)
for i in range(ndim):
    fit_params.append([pr[1,i], pr[2,i]-pr[1,i], pr[1,i]-pr[0,i],(pr[2,i]-pr[0,i])/2])

fit_params = numpy.array(fit_params)
print(fit_params)

fit_g = MCMC_program.shear_slope(fit_params[:,0], ra, dec)

h5f = h5py.File(parent_path + "result/%s_result_expo_%d_slope.hdf5"%(shear_cmd, expo_num),"w")
h5f["/chain"] = sampler.chain
h5f["/params"] = fit_params
h5f["/ra"] = ra
h5f["/dec"] = dec
h5f["/shear"] = fit_g
# h5f["/g1"] = g1
# h5f["/g2"] = g2
h5f.close()

inverse = range(nx-1, -1, -1)

my, mx = numpy.mgrid[-half_side:half_side,-half_side:half_side]
y, x = my*pixel_scale, mx*pixel_scale
shear_field = MCMC_program.shear_slope(parameters, x, y)[inverse]
shear_field_pdf = MCMC_program.shear_slope(fit_params[:, 0], x, y)[inverse]

g_min = min(shear.min(), fit_g.min(), shear_field_pdf.min())
g_max = max(shear.max(), fit_g.max(), shear_field_pdf.max())

img = Image_Plot()
img.subplots(2, 2)

cmap_g = plt.get_cmap('jet')
norm_g = plt.Normalize(vmin=g_min, vmax=g_max)

img.axs[0][0].scatter(ra, dec, color=cmap_g(norm_g(shear)), s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][0])

img.axs[0][1].scatter(ra, dec, color=cmap_g(norm_g(fit_g)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][1])

fig = img.axs[1][0].imshow(shear_field_pdf, vmin=g_min, vmax=g_max,cmap="jet")
img.figure.colorbar(fig, ax=img.axs[1][0])

fig = img.axs[1][1].imshow(shear_field-shear_field_pdf, cmap="jet")
img.figure.colorbar(fig, ax=img.axs[1][1])

yfmt = ["%.1f" % (8 - 4 * k) for k in range(5)]
xfmt = ["%.1f" % (-8 + 4 * k) for k in range(5)]

# img.set_ticklabel(1, 0, 0, len(xfmt), xfmt)
# img.set_ticklabel(1, 0, 1, len(xfmt), xfmt)
#
# img.set_ticklabel(1, 1, 0, len(xfmt), xfmt)
# img.set_ticklabel(1, 1, 1, len(xfmt), xfmt)

for i in range(2):
    for j in range(2):
        img.set_ticklabel(i, j, 0, len(yfmt), yfmt)
        img.set_ticklabel(i, j, 1, len(xfmt), xfmt)
        img.set_label(i, j, 0, "DEC.[arcmin]")
        img.set_label(i, j, 1, "R.A.[arcmin]")
img.save_img(parent_path + "pic/%s_mcmc_fit_expo_%d_slope.png"%(shear_cmd, expo_num))
img.close_img()




