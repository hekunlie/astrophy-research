from sys import path
import os
my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
path.append('%s/work/mylib/'%my_home)
path.append('D:/GitHub/astrophy-research/mylib')
path.append('D:/GitHub/astrophy-research/multi_shear_detect')
import numpy
import gauss_fit_fun
from plot_tool import Image_Plot
import dynesty
from dynesty import plotting as dyplot
import matplotlib.pyplot as plt



x = numpy.linspace(-10,10,501)

# # linear
# a1,a2 = 4.12, -1.55
# fx = a1*x+a2
#
# def logproir_linear(theta):
#     """log-likelihood"""
#     a1, a2 = theta
#     f1 = a1*x + a2
#     return -0.5*(numpy.sum((f1-fx)**2))
#
# def prior_transform_linear(utheta):
#     uw, umu = utheta
#     w = uw*20 - 10
#     mu = umu*10 - 5
#     return w, mu
#
# truths = [a1,a2]
# labels = ["$a_1$","$a_2$"]
#
# sampler = dynesty.NestedSampler(logproir_linear,prior_transform_linear,ndim=2,nlive=1000)
# sampler.run_nested()
# result = sampler.results
#
#
# fig, axes = dyplot.cornerplot(result, truths=truths, show_titles=True, title_kwargs={'y': 1.04}, labels=labels,
#                               fig=plt.subplots(2, 2, figsize=(15, 15)))
#
# plt.savefig("linear_1.png")
# plt.show()
#
#
# a1, a2, a3 = 4.12, -1.55, 5
# fx = a1*x**2 + a2*x + a3
#
# def logproir_linear(theta):
#     """log-likelihood"""
#     a1, a2, a3 = theta
#     f1 = a1*x**2 + a2*x + a3
#     return -0.5*(numpy.sum((f1-fx)**2))
#
# def prior_transform_linear(utheta):
#     ua1, ua2, ua3 = utheta
#     a1 = ua1*20 - 10
#     a2 = ua2*20 - 10
#     a3 = ua3*20 - 10
#     return a1,a2,a3
#
# truths = [a1,a2,a3]
# labels = ["$a_1$","$a_2$","$a_3$"]
#
# sampler = dynesty.NestedSampler(logproir_linear,prior_transform_linear,ndim=3,nlive=1000)
# sampler.run_nested()
# result = sampler.results
#
#
# fig, axes = dyplot.cornerplot(result, truths=truths, show_titles=True, title_kwargs={'y': 1.04}, labels=labels,
#                               fig=plt.subplots(3, 3, figsize=(15, 15)))
# plt.savefig("linear_2.png")
# plt.show()



# Gaussian

# mu,sig = 1, 2.3
# # single Gaussina test
# fx = gauss_fit_fun.gauss(x,mu,sig)
# def logproir_1(theta):
#     """log-likelihood"""
#     mu,sig = theta
#     f1 = gauss_fit_fun.gauss(x, mu,sig)
#     return -0.5*numpy.sum(1000*(f1-fx)**2)
#
# def prior_transform_1(utheta):
#     umu, usig = utheta
#     mu = umu*10 - 5
#     sig = usig*5
#     return mu, sig
#
# truths = [mu,sig]
# labels = ["$\mu$","$\sigma$"]
#
# sampler = dynesty.NestedSampler(logproir_1, prior_transform_1, ndim=2, nlive=500)
# sampler.run_nested(dlogz=0.01)
# result = sampler.results
#
# fig, axes = dyplot.cornerplot(result, truths=truths, show_titles=True, title_kwargs={'y': 1.04}, labels=labels,
#                               fig=plt.subplots(2, 2, figsize=(10, 10)))
# plt.savefig("gauss_1.png")
# plt.show()
#
# paras = result.samples
# para_1 = paras[:,0]
# para_2 = paras[:,1]
#
# loglikelihood = result.logl
# idx = loglikelihood > -0.1
#
#
# img = Image_Plot()
# img.subplots(1,3)
#
# img.axs[0][0].hist(loglikelihood[idx],100)
# img.axs[0][1].hist(para_1[idx],100)
# img.axs[0][2].hist(para_2[idx],100)
# plt.show()


# # test single Gaussina * line
# mu,sig = 1, 2.3
# a1,a2 = 0.31, -1.2
#
# fx1 = a1*x**2+a2
# fx2 = gauss_fit_fun.gauss(x,mu,sig)
# fx = fx1*fx2
#
# img = Image_Plot()
# img.subplots(1,1)
# # img.axs[0][0].plot(x,fx1)
# # img.axs[0][0].plot(x,fx2)
# img.axs[0][0].plot(x,fx)
# img.show_img()
#
#
# def logproir_1(theta):
#     """log-likelihood"""
#     mu,sig,a1,a2 = theta
#     f1 = gauss_fit_fun.gauss(x, mu,sig)*(a1*x**2+a2)
#     return -0.5*numpy.sum(1000*(f1-fx)**2)
#
# def prior_transform_1(utheta):
#     umu, usig,ua1,ua2 = utheta
#
#     mu = umu*10 - 5
#     sig = usig*5
#     a1 = 10*ua1 - 5
#     a2 = 10*ua2 - 5
#
#     return mu, sig, a1, a2
#
# truths = [mu,sig,a1,a2]
# labels = ["$\mu$","$\sigma$","$a_1$","$a_2$"]
#
# sampler = dynesty.NestedSampler(logproir_1, prior_transform_1, ndim=4, nlive=500)
# sampler.run_nested(dlogz=0.01)
# result = sampler.results
#
# fig, axes = dyplot.cornerplot(result, truths=truths, show_titles=True, title_kwargs={'y': 1.04}, labels=labels,
#                               fig=plt.subplots(4, 4, figsize=(16, 16)))
# plt.savefig("gauss_1_linear.png")
# plt.show()



# test two Gaussian
w1, mu1, sig1 = 5, 1.5, 4
w2, mu2, sig2 = 2, -2, 5

gauss_1 = gauss_fit_fun.gauss_coeff(x,w1,mu1,sig1)
gauss_2 = gauss_fit_fun.gauss_coeff(x,w2,mu2,sig2)
gauss_total = gauss_1 + gauss_2

def logproir_2(theta):
    """log-likelihood"""
    lw1, lmu1, lsig1, lw2, lmu2, lsig2= theta
    # lw1, lmu1, lsig1, lw2, lmu2, lsig2, lnf = theta
    f1 = gauss_fit_fun.gauss_coeff(x, lw1, lmu1, lsig1)
    f2 = gauss_fit_fun.gauss_coeff(x, lw2, lmu2, lsig2)
    model = f1 + f2
    return -0.5*(numpy.sum((model-gauss_total)**2))
    # inv_sig = 1./(model**2*numpy.exp(2*lnf))
    # return -0.5*(numpy.sum(((model-gauss_total)**2)/inv_sig - numpy.log(inv_sig)))

def prior_transform_2(utheta):
    uw1, umu1, usig1, uw2, umu2, usig2 = utheta
    pw1 = 10*uw1
    pmu1 = umu1*10 - 5
    psig1 = usig1*10
    pw2 = 10*uw2
    pmu2 = umu2*10 - 5
    psig2 = usig2*10
    # plnf = ulnf*20-10
    return pw1,pmu1, psig1,pw2,pmu2,psig2#,plnf

truths = [w1, mu1, sig1, w2, mu2, sig2]
labels = ["$w_1$", "$\mu_1$", "$\sigma_1$", "$w_1$","$\mu_2$", "$\sigma_2$"]

sampler = dynesty.NestedSampler(logproir_2, prior_transform_2, ndim=6, nlive=2000)
sampler.run_nested(dlogz=0.01)

# sampler = dynesty.DynamicNestedSampler(logproir_2, prior_transform_2, ndim=6, nlive=2000)
# sampler.run_nested(dlogz_init=0.05)
result = sampler.results

fig, axes = dyplot.cornerplot(result, truths=truths, show_titles=True, title_kwargs={'y': 1.04}, labels=labels,
                              fig=plt.subplots(6, 6, figsize=(15, 15)))
plt.savefig("gauss_2.png")
plt.show()
