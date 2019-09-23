import numpy
from sys import path,argv
path.append("/home/hklee/work/mylib/")
path.append("/home/hkli/work/mylib/")
path.append("/home/hklee/work/multi_shear_detect/")
import tool_box
from plot_tool import Image_Plot
from Fourier_Quad import Fourier_Quad
import gauss_fit_fun
import scipy
import copy
import h5py
from mpi4py import MPI


comm = MPI.COMM_WORLD
rank = comm.Get_rank()


fq = Fourier_Quad(12,1243)

# multi-signal
rng = numpy.random.RandomState(int(argv[1]) + int(rank*1234))
total_num = 20000000
signal_ini = [0.04, -0.04]
sigma = 0.06
fit_try = 10000
chisq_num = 100
gh = numpy.linspace(-0.1, 0.1, chisq_num)
dg = gh[1] - gh[0]

mix_data = numpy.zeros((total_num,))
for i in range(len(signal_ini)):
    num_sub = int(total_num/len(signal_ini))
    mix_data[i*num_sub:(i+1)*num_sub] = rng.normal(signal_ini[i], sigma, num_sub)

signals = copy.deepcopy(signal_ini)
bin_num = 8

# the signal of the mixture
ghat, ghat_sig = fq.find_shear(mix_data, 1, 8, left=-0.1, right=0.1)[:2]
signals.append(0)
print(signals, ghat, ghat_sig)


# calculate the \chi^2
# chisq = gauss_fit_fun.cal_chisq(mix_data, 1, bin_num, gh)

# calculate the flow
num_move = gauss_fit_fun.get_flow(mix_data, 1, gh)
text_content = []

print("One Gaussian fitting")
############################################################################
# fit
x_fit, y_fit = (gh[1:] + gh[:chisq_num-1])/2, num_move/num_move.sum()/dg

# fit one component
err_1 = 1.e20
err_1s = []
paras_1 = []
for i in range(fit_try):
    try:
        p0 = rng.uniform(0, 10, 3).tolist()
        fit_res_1 = scipy.optimize.curve_fit(gauss_fit_fun.gauss_coeff, x_fit, y_fit,
                                             p0=p0, bounds=([0, -10, 0], [5., 10, 10]))[0]
        err_ = gauss_fit_fun.fit_chisq_1_coeff(fit_res_1, x_fit, y_fit)
        if err_ <= err_1:
            a, b, c = fit_res_1
            err_1 = err_
        err_1s.append(err_)
        paras_1.append(fit_res_1)
    except:
        pass

err_1s = numpy.array(err_1s)
paras_1 = numpy.array(paras_1)


# ############################################################################
# fit two components

print("Two Gaussian's fitting")
err_2 = 1.e20
err_2s = []
paras_2 = []
for i in range(fit_try):
    try:
        p0 = rng.uniform(0, 10, 6).tolist()
        fit_res_2 = scipy.optimize.curve_fit(gauss_fit_fun.gauss_2_coeff, x_fit, y_fit,
                                             p0=p0, bounds=([0., -10, 0, 0., -10, 0], [5., 10, 10, 5., 10, 10]))[0]
        err_ = gauss_fit_fun.fit_chisq_2_coeff(fit_res_2, x_fit, y_fit)
        if err_ < err_2:
            num1, mu1, sig1, num2, mu2, sig2 = fit_res_2
            err_2 = err_
        err_2s.append(err_)
        paras_2.append(fit_res_2)
    except:
        pass
err_2s = numpy.array(err_2s)
paras_2 = numpy.array(paras_2)

h5f = h5py.File("pic/cache/cache_%d.hdf5"%rank, "w")
h5f["/chisq_1"] = err_1s
h5f["/para_1"] = paras_1
h5f["/chisq_2"] = err_2s
h5f["/para_2"] = paras_2
h5f["/signal"] = numpy.array(signal_ini)
h5f["/sigma"] = numpy.array([sigma])
h5f["/x_fit"] = x_fit
h5f["/y_fit"] = y_fit
h5f.close()





