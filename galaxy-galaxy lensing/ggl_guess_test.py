# import matplotlib
from sys import path,argv
path.append("/home/hklee/work/mylib")
from hk_plot_tool import Image_Plot
import hk_tool_box
import numpy
import h5py
import hk_FQlib
import time
# import hk_c4py
# import os
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
cpus = comm.Get_size()


def cal_signal(guess, data, pdf_bin):
    chisq = numpy.zeros_like(guess)

    bin_num = pdf_bin.shape[0] - 1
    bin_num2 = int(bin_num / 2)

    for i in range(guess.shape[0]):
        num = numpy.histogram(data - guess[i], pdf_bin)[0]
        #         print(num)
        n1, n2 = numpy.flip(num[:bin_num2], axis=0), num[bin_num2:]
        chisq[i] = numpy.sum((n1 - n2) ** 2 / (n1 + n2)) / 2

    idx = chisq <= chisq.min() + 150
    x, y = guess[idx], chisq[idx]
    #     img = Image_Plot()
    #     img.subplots(1,1)
    #     img.axs[0][0].plot(guess,chisq)
    #     img.show_img()
    coeff = tool_box.fit_1d(x, y, 2, "scipy")

    # y = a1 + a2*x + a3*x^2 = a3(x+a2/2/a3)^2 +...
    # gh = - a2/2/a3, gh_sig = 1/ sqrt(1/2/a3)
    gh = -coeff[1] / 2. / coeff[2]
    gh_sig = 0.70710678118 / numpy.sqrt(coeff[2])

    return gh, gh_sig, guess, chisq




guess_num = 200
sym_guess_m = 100
asym_guess_m = 50
sym_gh = numpy.zeros((guess_num, ))
asym_gh = numpy.zeros((guess_num, ))

sym_gh[:sym_guess_m] = -tool_box.set_bin_log(0.01, 200, sym_guess_m)
sym_gh[sym_guess_m:] = tool_box.set_bin_log(0.01, 200, sym_guess_m)
sym_gh = numpy.sort(sym_gh)

asym_gh[:asym_guess_m] = -tool_box.set_bin_log(0.01, 200, asym_guess_m)
asym_gh[asym_guess_m:] = tool_box.set_bin_log(0.01, 200, guess_num - asym_guess_m)
asym_gh = numpy.sort(asym_gh)

# 2 bins
bin_num = int(argv[1])
realization_num = 200
shear = [20, 140, 20, 140]
sigma = 1500
data_num = [20000, 20000, 1000000, 1000000]


rng = numpy.random.RandomState(12314)
rng_ch = numpy.random.RandomState(int(1231 + rank*100))

realization_list = [i for i in range(realization_num)]

my_list = tool_box.alloc(realization_list, cpus)[rank]


itemsize = MPI.DOUBLE.Get_size()
if rank == 0:
    # bytes for 10 double elements
    nbytes = int(2*len(shear)* realization_num)*itemsize
else:
    nbytes = 0

win1 = MPI.Win.Allocate_shared(nbytes, itemsize, comm=comm)
buf1, itemsize = win1.Shared_query(0)
results = numpy.ndarray(buffer=buf1, dtype='d', shape=(2*len(shear),realization_num))



for i in range(len(shear)):
    data = rng.normal(shear[i], sigma, data_num[i])
    for j in my_list:

        data_ch = rng_ch.choice(data, int(0.1 * data_num[i]), replace=False)

        pdf_bin = FQlib.set_bin(data_ch, bin_num, 1000000)

        gh = cal_signal(sym_gh, data, pdf_bin)[0]
        results[i, j] = gh

        gh = cal_signal(asym_gh, data, pdf_bin)[0]
        results[i+len(shear), j] = gh

comm.Barrier()
print(results)


if rank == 0:

    img = Image_Plot()
    img.subplots(2, 2)
    for i in range(2):
        for j in range(2):
            tag = int(i * 2 + j)
            st, ed = min(results[tag].min(), results[tag+len(shear)].min()), max(results[tag+len(shear)].max(), results[tag+len(shear)].max())
            hist_bin = numpy.linspace(st, ed, 50)

            img.axs[i][j].hist(results[tag], hist_bin, label="sym guess")
            img.axs[i][j].hist(results[tag + len(shear)], hist_bin, label="asym guess")

            img.axs[i][j].legend()
    img.save_img("/home/hklee/work/test_%dbin.pdf"%bin_num)
    img.close_img()
comm.Barrier()