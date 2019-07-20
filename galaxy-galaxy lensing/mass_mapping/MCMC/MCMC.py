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
    a1, a2, a3, a4, a5 = theta
    if -0.1 < a1 < 0.1 and -0.1 < a2 < 0.1 and -0.1 < a3 < 0.1 and -0.1 < a4 < 0.1 and -0.1 < a5 < 0.1:
        return 0.0
    return -numpy.inf

def ln_prob(theta, G, bins, bin_num2, inverse, x, x2, x3, x4, signal_num):
    lp = ln_gh_prior(theta)
    if not numpy.isfinite(lp):
        return -numpy.inf
    else:
        a1, a2, a3, a4, a5 = theta
        G_h = G - a1 - a2*x - a3*x2 - a4*x3 - a5*x4
        xi = 0
        for i in range(signal_num):
            num = numpy.histogram(G_h[i], bins)[0]
            n1 = num[0:bin_num2][inverse]
            n2 = num[bin_num2:]
            xi += numpy.sum((n1 - n2) ** 2 / (n1 + n2))*0.5
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


# a1 = -0.02
# a2 = -0.02
# a3 = 0.05

a1 = -0.035
a2 = 0.01
a3 = 0.02
a4 = 0
a5 = 0

num = int(argv[1])
ncpus = int(argv[2])
signal_num = 15

nwalkers = 300
ndim = 5
step = 600
print("Walker: %d. Step: %d."%(nwalkers, step))
fq = Fourier_Quad(10, 112)
bin_num = 8
bin_num2 = int(bin_num / 2)
x = numpy.linspace(-1, 1, signal_num).reshape((signal_num, 1))
x2 = x*x
x3 = x*x*x
x4 = x*x*x*x
signals = a1 + a2*x + a3*x2 + a4*x3 + a5*x4
parameters = [a1, a2, a3, a4, a5]
print("Signals: ", signals[:,0],".\n%.4f + %.4f*x + %.4f*x^2 + %.4f*x^3 + %.4f*x^4"%(a1,a2,a3, a4, a5))

ellip = numpy.zeros((signal_num, num))

img = Image_Plot()
img.subplots(1,1)

fq_shear = numpy.zeros((2,signal_num))
for i in range(signal_num):
    # rng = numpy.random.RandomState(i+1)
    # ellip[i] = rng.normal(signals[i,0], 0.3, num)
    ellip[i] = numpy.random.normal(signals[i,0], 0.3, num)
    # noise = rng.normal(0, numpy.abs(ellip[i])/5)
    # ellip[i] += noise
    t1 = time.time()
    gh, gh_sig = fq.fmin_g_new(ellip[i], numpy.ones_like(ellip[i]), 8)[:2]
    fq_shear[0, i] = gh
    fq_shear[1, i] = gh_sig
    t2 = time.time()
    print("signal:[%.4f at %.4f] %.4f (%.4f) [%d gal], Time: %.2f sec"%(signals[i,0], x[i,0], gh, gh_sig, num, t2-t1))
    img.axs[0][0].hist(ellip[i], 100, histtype="step", label="%.4f" % signals[i])
img.save_img("./pic/data_hist.png")
# img.show_img()
img.close_img()

# find the signal
all_ellip = ellip.reshape((num*signal_num,))
ellip_bins = fq.set_bin(all_ellip, 8, 1.5)
inverse = range(bin_num2 - 1, -1, -1)
print("Bins:", ellip_bins)


p0 = numpy.random.uniform(-0.02, 0.02, ndim*nwalkers).reshape((nwalkers, ndim))

t1 = time.time()
with Pool(ncpus) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, ln_prob, args=[ellip, ellip_bins, bin_num2, inverse, x, x2, x3, x4, signal_num],pool=pool)
    t2 = time.time()
    pos, prob, state = sampler.run_mcmc(p0, step)
t3 = time.time()
print("Time: %.2f sec, %2.f sec"%(t2-t1, t3-t2))


img = Image_Plot(fig_x=16, fig_y=4)
img.subplots(ndim,1)
for i in range(nwalkers):
    for j in range(ndim):
        img.axs[j][0].plot(range(step),sampler.chain[i, :, j], color='grey',alpha=0.6)
        img.axs[j][0].plot([0,step], [parameters[j], parameters[j]])
img.save_img("./pic/mcmc_walkers_nw_%d_stp_%d.png"%(nwalkers, step))
img.close_img()

samples = sampler.chain[:, 150:, :].reshape((-1, ndim))
print(samples.shape)
corner_fig = plt.figure(figsize=(10, 10))
fig = corner.corner(samples, labels=["$a_1$", "$a_2$", "$a_3$", "$a_4$", "$a_5$"], truths=[a1, a2, a3, a4, a5],
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".4f", title_kwargs={"fontsize": 12})
fig.savefig("./pic/mcmc_panel_nw_%d_stp_%d.png"%(nwalkers, step))

fit_params = []
pr = numpy.percentile(samples, [16, 50, 84], axis=0)
for i in range(ndim):
    fit_params.append([pr[1,i], pr[2,i]-pr[1,i], pr[1,i]-pr[0,i]])
fit_params_ = map(lambda v: (v[1], v[2] - v[1], v[1] - v[0]),zip(*numpy.percentile(samples, [16, 50, 84], axis=0)))

print(fit_params)
for para in fit_params_:
    print("%8.5f [%8.5f, %8.5f]"%(para[0], para[1], para[2]))


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



# for i in range(ndim):
#     img = Image_Plot()
#     img.subplots(1, 1)
#     img.axs[0][0].hist(sampler.flatchain[:, 0], 100, histtype="step", color='k')
#     img.save_img("mcmc_chisq.png")
#     # img.show_img()
#     img.close_img()
# pool.close()