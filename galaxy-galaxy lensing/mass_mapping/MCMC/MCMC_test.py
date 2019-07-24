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

nwalkers = 200
ndim = 4
step = 1000
print("Walker: %d. Step: %d."%(nwalkers, step))
fq = Fourier_Quad(10, 112)
bin_num = 8
bin_num2 = int(bin_num / 2)

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
shear = data[:, 5]
g1 = data[:, 6]
g2 = data[:, 7]
ra = data[:, 8]
dec = data[:, 9]

radius = numpy.sqrt(ra**2 + dec**2)
shear_n = shear + numpy.random.normal(0, 0.0005, len(shear))

sin_theta = ra / radius
cos_theta = dec / radius
sin_2theta = 2 * sin_theta * cos_theta
cos_2theta = cos_theta ** 2 - sin_theta ** 2

mg_bins = fq.set_bin(mg1, bin_num, 5)
inverse = range(bin_num2 - 1, -1, -1)
print("Bins:", mg_bins)

p0 = numpy.zeros((nwalkers, ndim))
p0[:, 0] = numpy.random.uniform(0, 1, nwalkers)
p0[:, 1] = numpy.random.uniform(-0.1, 0.1, nwalkers)
p0[:, 2] = numpy.random.uniform(-0.1, 0.1, nwalkers)
p0[:, 3] = numpy.random.uniform(1, 5, nwalkers)


print("Begin the test")

# MCMC
t1 = time.time()
with Pool(ncpus) as pool:
    sampler_test = emcee.EnsembleSampler(nwalkers, ndim, MCMC_program.ln_prob_g_test, args=[shear_n, ra, dec, 0], pool=pool)
    pos, prob, state = sampler_test.run_mcmc(p0, step)
samples = sampler_test.chain[:, 150:, :].reshape((-1, ndim))
print(samples.shape)

# Plot the walkers
img = Image_Plot(fig_x=16, fig_y=4)
img.subplots(ndim,1)
for i in range(nwalkers):
    for j in range(ndim):
        img.axs[j][0].plot(range(step),sampler_test.chain[i, :, j], color='grey',alpha=0.6)
        img.axs[j][0].plot([0,step], [parameters[j], parameters[j]])
img.save_img(parent_path + "pic/test_mcmc_walkers_nw_%d_stp_%d_expo_%d.png"%(nwalkers, step, expo_num))
img.close_img()

# Plot the contour
corner_fig = plt.figure(figsize=(10, 10))
fig = corner.corner(samples, labels=["$A$", "$\delta RA$", "$\delta DEC$", "$\sigma$"], truths=parameters,
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".4f", title_kwargs={"fontsize": 12})
fig.savefig(parent_path + "pic/test_mcmc_panel_nw_%d_stp_%d_expo_%d.png"%(nwalkers, step,expo_num))

# MCMC results
fit_params_mcmc = []
pr = numpy.percentile(samples, [16, 50, 84], axis=0)
for i in range(ndim):
    fit_params_mcmc.append([pr[1,i], pr[2,i]-pr[1,i], pr[1,i]-pr[0,i]])
fit_params_mcmc = numpy.array(fit_params_mcmc)
fit_g_mcmc = MCMC_program.shear_profile(fit_params_mcmc[:,0], ra, dec)[0]
t2 = time.time()
print(t2-t1)
print("MCMC:", fit_params_mcmc[:,0])

# least square results
t1 = time.time()
p0 = numpy.array([2, 0.1, 5, 10])
xy = numpy.row_stack((ra, dec))
fit_param_lst = least_squares(MCMC_program.surface_fit, p0, args=[xy, shear_n]).x
fit_g_lst = MCMC_program.shear_profile(fit_param_lst, ra, dec)[0]
t2 = time.time()
print("Least Sqaure: ", fit_param_lst)
print(t2-t1)


# residula
residual_mcmc = fit_g_mcmc - shear
residual_lst = fit_g_lst - shear

# Plot the fitted surface
g_min = min([shear.min(), fit_g_mcmc.min(), fit_g_lst.min()])
g_max = max([shear.max(), fit_g_mcmc.max(), fit_g_lst.max()])

diff_min = min(residual_lst.min(), residual_mcmc.min())
diff_max = max(residual_lst.max(), residual_mcmc.max())

img = Image_Plot()
img.subplots(2, 3)

cmap_g = plt.get_cmap('jet')
norm_g = plt.Normalize(vmin=g_min, vmax=g_max)

img.axs[0][0].set_title("Input")
img.axs[0][0].scatter(ra, dec, color=cmap_g(norm_g(shear)), s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][0])

img.axs[0][1].set_title("MCMC")
img.axs[0][1].scatter(ra, dec, color=cmap_g(norm_g(fit_g_mcmc)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][1])

img.axs[0][2].set_title("Least Square")
img.axs[0][2].scatter(ra, dec, color=cmap_g(norm_g(fit_g_lst)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[0][2])


# residual
cmap_g = plt.get_cmap('jet')
norm_g = plt.Normalize(vmin=diff_min, vmax=diff_max)

img.axs[1][1].set_title("$\Delta$ MCMC")
img.axs[1][1].scatter(ra, dec, color=cmap_g(norm_g(residual_mcmc)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[1][1])

img.axs[1][2].set_title("$\Delta$ Least square")
img.axs[1][2].scatter(ra, dec, color=cmap_g(norm_g(residual_lst)),s=1)
sm = plt.cm.ScalarMappable(cmap=cmap_g, norm=norm_g)
sm._A = []
plt.colorbar(sm, ax=img.axs[1][2])

img.save_img(parent_path + "pic/test_mcmc_fit_nw_%d_stp_%d_expo_%d.png"%(nwalkers, step, expo_num))
img.close_img()



