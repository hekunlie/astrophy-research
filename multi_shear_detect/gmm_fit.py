import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy
from sys import path
path.append("E:/Github/astrophy-research/mylib")
path.append("/home/hklee/work/mylib")
from plot_tool import Image_Plot
import scipy
from sklearn.mixture import GaussianMixture
import time
import emcee
import corner
from multiprocessing import Pool


def gauss(params, x):
    a2, a3 = params
    if a2 < 0:
        print("Negative value",a2)
    f1 = 1/numpy.sqrt(2*numpy.pi*(a2**2))*numpy.exp(-(x-a3)**2/2/a2**2)
    # a2, a3 = p
    # f1 = 1./numpy.sqrt(2*numpy.pi)/a2*numpy.exp(-(x-a3)**2/2/a2**2)
    dx = x[1] - x[0]
    return f1, dx


def mix_gauss(params, x):
    w1, w2 = params[0], params[3]
    f1 = gauss(params[1:3], x)[0]
    f2 = gauss(params[4:], x)[0]
    dx = x[1] - x[0]
    denominator = numpy.sum(w1*f1*dx+w2*f2*dx)
    if denominator == 0:
        print("Divided by ZERO")
    return (w1*f1 + w2*f2)/denominator, f1, f2

def error(p, x, y):
    return mix_gauss(p, x)[0] - y

def para_prior(params):
    tag = 0
    num = len(params)
    fit_range = [[0, 1], [0, 30], [-30, 30], [0, 1], [0, 30], [-30, 30]]
    for i in range(num):
        if fit_range[i][0] < params[i] < fit_range[i][1]:
            tag += 1
    if tag == num and abs(params[0] + params[3] - 1) <= 0.00001:
        return 0
    else:
        return -numpy.inf

def chi_sq(params, x, y):
    lp = para_prior(params)
    if not numpy.isfinite(lp):
        return -numpy.inf
    else:
        return lp-numpy.sum((mix_gauss(params, x)[0] - y)**2*10000)


def str_format(para_list, sub_len):
    strs = ""
    n = int(len(para_list)/sub_len)
    for i in range(n):
        for j in range(sub_len):
            strs += "%6.3f      "%para_list[i*sub_len + j]
    return strs


# SAMPLE
seed = 1279734
total_num = 600000
num1 = 300000
rng = numpy.random.RandomState(seed)
# sig, mean, sig, mean
sig1, mean1 = 4, -5
sig2, mean2 = 4, 2
w1_true, w2_true = num1*1.0/total_num, (total_num - num1*1.0)/total_num
para_true = [w1_true, sig1, mean1, w2_true, sig2, mean2]
print("The input:")
print("weight    sigma     mean    weight      sigma       mean")
print(str_format(para_true,3), "\n")
ini_samples = numpy.zeros((total_num,))
ini_samples[:num1] = rng.normal(mean1, sig1, num1)
ini_samples[num1:] = rng.normal(mean2, sig2, total_num-num1)
# idx1 = samples >= -20
# idx2 = samples <= 20
# samples = samples[idx1&idx2]

# histogram
bin_num = 200
fx, bins = numpy.histogram(ini_samples, bin_num, normed=True)[:2]
x = (bins[:bin_num] + bins[1:])/2
dx = x[1] - x[0]

f_mix, f1_true, f2_true = mix_gauss(para_true, x)

denominator = numpy.sum(w1_true*f1_true*dx + w2_true*f2_true*dx)
f1 = w1_true*f1_true/denominator
f2 = w2_true*f2_true/denominator
# print(numpy.sum(f_mix*dx))
# funs = [f_mix, f1, f2, f1+f2]
# img = Image_Plot()
# img.subplots(1,1)
# for i in range(4):
#     img.axs[0][0].plot(x, funs[i])
# img.show_img()
# img.close()
# print(numpy.sum(fx*dx))
# print(fx.shape, x.shape)

# #####################################   Least squares   ################################################################
# # Least squares
# coeff_0 = [1.1, 1.1, 0.1, 5, 1.5, 1.1]
# res = scipy.optimize.least_squares(error, coeff_0, args=(x, fx))
# # res = scipy.optimize.minimize(chi_sq, coeff_0, args=(x, fx))
# res_fit = res.x
#
# dn = numpy.abs(res_fit[0] + res_fit[3])
# w1, w2 = numpy.abs(res_fit[0])/dn, numpy.abs(res_fit[3])/dn
# # w1, w2 = numpy.abs(res_fit[0]), numpy.abs(res_fit[1])
#
# para_fit = [w1,res_fit[1],res_fit[2], w2, res_fit[4], res_fit[5]]
# f_fit, f1_fit_, f2_fit_ = mix_gauss(para_fit, x)
# denominator = numpy.sum(w1*f1_fit_*dx + w2*f2_fit_*dx)
# f1_fit = w1*f1_fit_/denominator
# f2_fit = w2*f2_fit_/denominator
# print("The leastsq:")
# print(res_fit)
# print("weight    sigma     mean    weight      sigma       mean")
# print(str_format(para_fit, 3), "\n")
# # print(numpy.sum(f_fit*dx))
# #####################################   Least squares   ################################################################

#####################################   GMM   ################################################################
gmm = GaussianMixture(n_components=2)
gmm.fit(ini_samples.reshape((-1,1)))
comp_label = gmm.predict(ini_samples.reshape((-1,1)))
# img = Image_Plot()
# img.subplots(1,1)
# img.axs[0][0].hist(comp_label,10)
# img.show_img()

idx = comp_label < 0.5
num1 = idx.sum()
idx = comp_label > 0.5
num2 = idx.sum()
w1 = num1/(num1+num2)
w2 = num2/(num1+num2)
# print(num1, num2)
means = gmm.means_[:,0]
covs = gmm.covariances_[:,0,0]
para_gmm = [w1, numpy.sqrt(covs)[0], means[0], w2, numpy.sqrt(covs)[1], means[1]]

f1_gmm_ = gauss(para_gmm[1:3], x)[0]
f2_gmm_ = gauss(para_gmm[4:], x)[0]

f_gmm = f1_gmm_*w1 + f2_gmm_*w2
denominator = numpy.sum(w1*f2_gmm_*dx + w2*f2_gmm_*dx)
f1_gmm = w1*f1_gmm_/denominator
f2_gmm = w2*f2_gmm_/denominator
print("The gmm:")
print("weight    sigma     mean    weight      sigma       mean")
print(str_format(para_gmm,3),"\n")
#####################################   GMM   ################################################################


#####################################   MCMC   ################################################################
ncpus = 10
nwalkers = 500
ndim = 6
step = 5000

p0 = numpy.zeros((nwalkers, ndim))
p0[:, 0] = numpy.random.uniform(0.0001, 1, nwalkers)
p0[:, 1] = numpy.random.uniform(0, 2, nwalkers)
p0[:, 2] = numpy.random.uniform(-10, 10, nwalkers)
p0[:, 3] = numpy.random.uniform(0.0001, 1, nwalkers)
p0[:, 4] = numpy.random.uniform(0, 2, nwalkers)
p0[:, 5] = numpy.random.uniform(-10, 10, nwalkers)

t1 = time.time()
with Pool(ncpus) as pool:
    sampler = emcee.EnsembleSampler(nwalkers, ndim, chi_sq, args=[x, fx], pool=pool)
    pos, prob, state = sampler.run_mcmc(p0, step)
samples = sampler.chain[:, 150:, :].reshape((-1, ndim))
print(samples.shape)

# Plot the walkers
img = Image_Plot(fig_x=18, fig_y=3)
img.subplots(ndim,1)
for i in range(nwalkers):
    for j in range(ndim):
        img.axs[j][0].plot(range(step),sampler.chain[i, :, j], color='grey',alpha=0.6)
        img.axs[j][0].plot([0,step], [para_true[j], para_true[j]])
img.save_img("mcmc_walker.png")
img.close_img()

# Plot the contour
corner_fig = plt.figure(figsize=(10, 10))
fig = corner.corner(samples, labels=["$w_1$", "$\sigma_1$", "$\mu_1$", "$w_2$", "$\sigma_2$", "$\mu_2$",], truths=para_true,
                    quantiles=[0.16, 0.5, 0.84], show_titles=True, title_fmt=".4f", title_kwargs={"fontsize": 12})
fig.savefig("mcmc_contour.png")

# MCMC results
params_mcmc = []
pr = numpy.percentile(samples, [16, 50, 84], axis=0)
for i in range(ndim):
    params_mcmc.append(pr[1,i])

f_mc, f1_mc_, f2_mc_ = mix_gauss(params_mcmc, x)
w1, w2 = params_mcmc[0], params_mcmc[3]
denominator = numpy.sum(w1*f1_mc_*dx + w2*f2_mc_*dx)
f1_mc = w1*f1_mc_/denominator
f2_mc = w2*f2_mc_/denominator

t2 = time.time()
print("The MCMC:")
print("weight    sigma     mean    weight      sigma       mean")
print(str_format(params_mcmc, 3),"\n")
print(t2-t1)
#####################################   MCMC   ################################################################



# Plot the results
lw = 2.5
img = Image_Plot()
img.set_style()
img.subplots(1,2)
img.axs[0][0].fill_between(x,0, fx, alpha=0.5, label="PDF",color="C9")
img.axs[0][0].plot(x,f_mix,label="True total", linewidth=lw)
img.axs[0][0].plot(x,f1, label="True f1", linestyle="dotted", linewidth=lw)
img.axs[0][0].plot(x,f2, label="True f2",linestyle="dotted", linewidth=lw)
img.axs[0][0].plot(x,f_gmm, label="GMM: total", linewidth=lw)
img.axs[0][0].plot(x,f1_gmm,label="GMM: f1",linestyle="--", linewidth=lw)
img.axs[0][0].plot(x,f2_gmm,label="GMM: f2",linestyle="--", linewidth=lw)

img.axs[0][1].fill_between(x,0, fx, alpha=0.5, label="PDF",color="C9")
img.axs[0][1].plot(x,f_mix,label="True total", linewidth=lw)
img.axs[0][1].plot(x,f1, label="True f1", linestyle="dotted", linewidth=lw)
img.axs[0][1].plot(x,f2, label="True f2",linestyle="dotted", linewidth=lw)

img.axs[0][1].plot(x,f_mc, label="MCMC: total", linewidth=lw)
img.axs[0][1].plot(x,f1_mc,label="MCMC: f1",linestyle="--", linewidth=lw)
img.axs[0][1].plot(x,f2_mc,label="MCMC: f2",linestyle="--", linewidth=lw)

# img.axs[0][1].plot(x,f_fit, label="Least_sq: total", linewidth=lw)
# img.axs[0][1].plot(x,f1_fit,label="Least_sq: f1",linestyle="--", linewidth=lw)
# img.axs[0][1].plot(x,f2_fit,label="Least_sq: f2",linestyle="--", linewidth=lw)

ys = img.axs[0][0].set_ylim()
dy = (ys[1] - ys[0])*0.3
img.axs[0][0].set_ylim((ys[0], ys[1]+dy))
ys = img.axs[0][1].set_ylim()
dy = (ys[1] - ys[0])*0.3
img.axs[0][1].set_ylim((ys[0], ys[1]+dy))

img.axs[0][0].legend(ncol=2, fontsize=img.legend_size-4,loc="upper left")
img.axs[0][1].legend(ncol=2, fontsize=img.legend_size-4,loc="upper left")
for i in range(2):
    img.set_label(0,i,0,"$P(X)$")
    img.set_label(0,i,1,"$X$")
img.save_img("fit.png")
# img.show_img()