from sys import path
import os
# my_home = os.popen("echo $MYWORK_DIR").readlines()[0][:-1]
# path.append('%s/work/mylib/'%my_home)
path.append('D:/GitHub/astrophy-research/mylib')
path.append('D:/GitHub/astrophy-research/multi_shear_detect')
import numpy
import gauss_fit_fun
from plot_tool import Image_Plot
import dynesty
from dynesty import plotting as dyplot
import matplotlib.pyplot as plt



x = numpy.linspace(-10,10,501)

# linear
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
# plt.show()

a1, a2, a3 = 4.12, -1.55, 5
fx = a1*x**2 + a2*x + a3

def logproir_linear(theta):
    """log-likelihood"""
    a1, a2, a3 = theta
    f1 = a1*x**2 + a2*x + a3
    return -0.5*(numpy.sum((f1-fx)**2))

def prior_transform_linear(utheta):
    ua1, ua2, ua3 = utheta
    a1 = ua1*20 - 10
    a2 = ua2*20 - 10
    a3 = ua3*20 - 10
    return a1,a2,a3

truths = [a1,a2,a3]
labels = ["$a_1$","$a_2$","$a_3$"]

sampler = dynesty.NestedSampler(logproir_linear,prior_transform_linear,ndim=3,nlive=1000)
sampler.run_nested()
result = sampler.results


fig, axes = dyplot.cornerplot(result, truths=truths, show_titles=True, title_kwargs={'y': 1.04}, labels=labels,
                              fig=plt.subplots(3, 3, figsize=(15, 15)))
plt.show()


exit()
# Gaussian
w1, mu1, sig1 = 0.4, 1.5, 4
w2, mu2, sig2 = 0.6, -2, 5

gauss_1 = gauss_fit_fun.gauss_coeff(x,w1,mu1,sig1)
gauss_2 = gauss_fit_fun.gauss_coeff(x,w2,mu2,sig2)
gauss_total = gauss_1 + gauss_2

# single Gaussina test
def logproir_1(theta):
    """log-likelihood"""
    w1, mu1, sig1 = theta
    f1 = gauss_fit_fun.gauss_coeff(x, w1, mu1, sig1)
    return -0.5*(numpy.sum((f1-gauss_1)**2))

def prior_transform_1(utheta):
    uw, umu, usig = utheta
    w = uw*2 - 1
    mu = umu*10 - 5
    sig = usig*5
    return w, mu, sig

truths = [w1,mu1,sig1]
labels = ["$w$","$\mu$","$\sigma$"]

sampler = dynesty.NestedSampler(logproir_1,prior_transform_1,ndim=3,nlive=1000)
sampler.run_nested()
result = sampler.results


fig, axes = dyplot.cornerplot(result, truths=truths, show_titles=True, title_kwargs={'y': 1.04}, labels=labels,
                              fig=plt.subplots(3, 3, figsize=(10, 10)))
plt.show()


# test two Gaussian
def logproir_2(theta,x,y):
    """log-likelihood"""
    w1,mu1, sig1, w2, mu2, sig2 = theta
    f1 = gauss_fit_fun.gauss_coeff(x, w1, mu1, sig1)
    f2 = gauss_fit_fun.gauss_coeff(x, w2, mu2, sig2)
    f = f1 + f2
    return -0.5*(numpy.sum((f-y)**2))

def prior_transform_2(theta):
    pass