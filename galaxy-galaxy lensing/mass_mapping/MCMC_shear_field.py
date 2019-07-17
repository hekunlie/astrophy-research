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


def ln_gh_prior(theta):
    tag = 0
    num = 4
    fit_range = [[0, 2],[-10, 10],[-10, 10],[0, 36]]
    for i in range(num):
        if fit_range[i][0] <= theta[i] <= fit_range[i][1]:
            tag += 1
    if tag == num:
        return 0.0
    else:
        return -numpy.inf


def ln_prob_g1(theta, G1, NU1, bins, bin_num2, inverse, x, y, cos_2theta, r):
    lp = ln_gh_prior(theta)
    if not numpy.isfinite(lp):
        return -numpy.inf
    else:
        g = theta[0] * r * numpy.exp(-((y - theta[1])**2 + (x - theta[2])**2)/2/theta[3])/numpy.pi/2/theta[3]
        G_h = G1 - NU1 * g * cos_2theta
        num = numpy.histogram(G_h[i], bins)[0]
        n1 = num[0:bin_num2][inverse]
        n2 = num[bin_num2:]
        xi = numpy.sum((n1 - n2) ** 2 / (n1 + n2))*0.5
        return lp - xi


def ln_prob_g2(theta, G2, NU2, bins, bin_num2, inverse, x, y, sin_2theta, r):
    lp = ln_gh_prior(theta)
    if not numpy.isfinite(lp):
        return -numpy.inf
    else:
        g = theta[0] * r * numpy.exp(-((y - theta[1])**2 + (x - theta[2])**2)/2/theta[3])/numpy.pi/2/theta[3]
        G_h = G2 + NU2 * g*sin_2theta
        num = numpy.histogram(G_h[i], bins)[0]
        n1 = num[0:bin_num2][inverse]
        n2 = num[bin_num2:]
        xi = numpy.sum((n1 - n2) ** 2 / (n1 + n2))*0.5
        return lp - xi


def result_fun(params, coord):
    f, f_sig = 0, 0
    tag = 0
    for para in params:
        x = coord**tag
        f += para[0]*x
        f_sig += (para[1] + para[2])/2*x
        tag += 1
    return f, f_sig


ncpus = int(argv[2])
expo_num = int(argv[1])
shear_cmd = argv[3]

nwalkers = 300
ndim = 5
step = 600
print("Walker: %d. Step: %d."%(nwalkers, step))
fq = Fourier_Quad(10, 112)
bin_num = 8
bin_num2 = int(bin_num / 2)

# number of grid
nx = 50
ny = nx
# the field size in unit of arcmin
delta_ra = 20
delta_dec = delta_ra
half_side = nx/2
# arcmin/pix
pixel_scale = delta_ra/nx
# parameters
sigma = 1.5 # arcmin
amplitude = 0.8
shift = (0, 0)
parameters = [amplitude, shift[0], shift[1], sigma**2]


parent_path = "/mnt/perc/hklee/CFHT/multi_shear/cluster_field/result/"
for i in range(expo_num):
    h5f = h5py.File(parent_path + "expo_%d.hdf5"%i,"r")
    temp = h5f["/data"].value
    h5f.close()
    if i == 0:
        data = temp
    else:
        data = numpy.row_stack((data, temp))
mg1 = data[:, 0]
mg2 = data[:, 1]
mnu1 = data[:, 2] + data[:, 3]
mnu2 = data[:, 2] - data[:, 3]
g1 = data[:, 5]
g2 = data[:, 6]
ra = data[:, 7]
dec = data[:, 8]
radius = numpy.sqrt(ra**2 + dec**2)
sin_theta = ra / radius
cos_theta = dec / radius
sin_2theta = 2 * sin_theta * cos_theta
cos_2theta = cos_theta ** 2 - sin_theta ** 2

mg_bins = fq.set_bin(mg1, 8, 5)
inverse = range(bin_num2 - 1, -1, -1)
print("Bins:", mg_bins)


p0 = numpy.random.uniform(-0.02, 0.02, ndim*nwalkers).reshape((nwalkers, ndim))

if shear_cmd == "g1":
    t1 = time.time()
    with Pool(ncpus) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_prob_g1,
                                        args=[mg1, mnu1, mg_bins, bin_num2, inverse, ra, dec, cos_2theta, radius],pool=pool)
        t2 = time.time()
        pos, prob, state = sampler.run_mcmc(p0, step)
    t3 = time.time()

else:
    t1 = time.time()
    with Pool(ncpus) as pool:
        sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_prob_g2,
                                        args=[mg2, mnu2, mg_bins, bin_num2, inverse, ra, dec, sin_2theta, radius],pool=pool)
        t2 = time.time()
        pos, prob, state = sampler.run_mcmc(p0, step)
    t3 = time.time()
print("Time: %.2f sec, %2.f sec" % (t2 - t1, t3 - t2))


img = Image_Plot(fig_x=16, fig_y=4)
img.subplots(ndim,1)
for i in range(nwalkers):
    for j in range(ndim):
        img.axs[j][0].plot(range(step),sampler.chain[i, :, j], color='grey',alpha=0.6)
        img.axs[j][0].plot([0,step], [parameters[j], parameters[j]])
img.save_img("./pic/%s_mcmc_walkers_nw_%d_stp_%d.png"%(shear_cmd, nwalkers, step))
img.close_img()

samples = sampler.chain[:, 150:, :].reshape((-1, ndim))
print(samples.shape)
corner_fig = plt.figure(figsize=(10, 10))
fig = corner.corner(samples, labels=["$A$", "$\delta RA$", "$\delta DEC$", "$\sigma^2$"], truths=parameters,
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".4f", title_kwargs={"fontsize": 12})
fig.savefig("./pic/%s_mcmc_panel_nw_%d_stp_%d.png"%(shear_cmd,nwalkers, step))

fit_params = []
pr = numpy.percentile(samples, [16, 50, 84], axis=0)
for i in range(ndim):
    fit_params.append([pr[1,i], pr[2,i]-pr[1,i], pr[1,i]-pr[0,i]])
fit_params_ = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),zip(*numpy.percentile(samples, [16, 50, 84], axis=0)))

print(fit_params)

mcmc_shear, mcmc_sig = result_fun(fit_params, x)

result = numpy.zeros((6, signal_num))
result[0] = x[:,0]
result[1] = signals[:,0]
result[2] = fq_shear[0]
result[3] = fq_shear[1]
result[4] = mcmc_shear[:,0]
result[5] = mcmc_sig[:,0]

h5f = h5py.File("./result.hdf5","w")
h5f["/chain"] = sampler.chain
h5f["/result"] = result
h5f.close()

img = Image_Plot()
img.subplots(1,2)

img.axs[0][0].plot(x, signals, color='k', label="True")
img.axs[0][0].errorbar(x, mcmc_shear,mcmc_sig, label="MCMC Recovered")
img.axs[0][0].errorbar(x, fq_shear[0], fq_shear[1], label="FQ Recovered")
img.set_label(0,0,0, "g")
img.set_label(0,0,1, "X")
img.axs[0][0].legend()

img.axs[0][1].plot(x, 100*(signals[:,0] - mcmc_shear[:,0]), label="MCMC: True - Recovered")
img.axs[0][1].plot(x, 100*(signals[:,0] - fq_shear[0]), label="FQ: True - Recovered")
img.set_label(0,1,0, "$10^2 \\times\Delta g$")
img.set_label(0,1,1, "X")
img.axs[0][1].legend()
img.subimg_adjust(0, 0.25)

img.save_img("./pic/mcmc_recover_nw_%d_stp_%d.png"%(nwalkers, step))
img.close_img()
