import numpy
from sys import path
path.append("/home/hklee/work/mylib/")
path.append("D:/GitHub/astrophy-research/mylib/")
from Fourier_Quad import Fourier_Quad



def gauss(x, mu, sig):
    """ gaussian PDF"""
    return 1./numpy.sqrt(2*numpy.pi*sig**2)*numpy.exp(-(x-mu)**2/2/sig/sig)

def gauss_coeff(x, coeff, mu, sig):
    """ gaussian PDF x coefficient, used for fitting the unnormalized histogram"""
    return coeff/numpy.sqrt(2*numpy.pi*sig**2)*numpy.exp(-(x-mu)**2/2/sig/sig)


def gauss_2(x, w1, mu1, sig1, mu2,sig2):
    """ two weighted gaussian PDFs"""
    f1 = w1*gauss(x, mu1, sig1)
    f2 = (1-w1)*gauss(x, mu2, sig2)
    return f1 + f2

def gauss_2_coeff(x, coeff1, mu1, sig1, coeff2, mu2, sig2):
    """ two weighted gaussian PDF x coefficient, used for fitting the unnormalized histogram"""
    f1 = gauss_coeff(x, coeff1, mu1, sig1)
    f2 = gauss_coeff(x, coeff2, mu2, sig2)
    return f1 + f2


def gauss_fun_2_(x, coeff1, coeff2):
    f1 = gauss(x, 0.04, 0.02)*coeff1
    f2 = gauss(x, -0.04, 0.02)*coeff2
    return f1 + f2


def fit_chisq_1(para, x, y):
    """ chi squared """
    return numpy.sum((gauss(x, para[0], para[1])-y)**2)

def fit_chisq_1_coeff(para,x,y):
    """ chi squared """
    return numpy.sum((gauss_coeff(x, para[0], para[1], para[2])-y)**2)


def fit_chisq_2(para, x, y):
    """ chi squared """
    a, b, c, d, e = para
    return numpy.sum((gauss_2(x, a, b, c, d, e)-y)**2)

def fit_chisq_2_coeff(para,x,y):
    """ chi squared """
    a, b, c, d, e, f = para
    return numpy.sum((gauss_2_coeff(x, a, b, c, d, e, f)-y)**2)


def get_flow(mg, mnu, gh):
    """ calculate the number of particle cross zero from right to left """
    gh_num = gh.shape[0]
    num_move = numpy.zeros((gh_num-1,))
    for i in range(gh_num-1):
        idx1 = mg - gh[i]*mnu >= 0
        idx2 = mg - gh[i+1]*mnu >= 0
        num_move[i] = idx1.sum() - idx2.sum()
    return num_move


def cal_chisq(mg, mnu, bin_num, gh):
    """ calculate \chi^2 of PDF_SYM """
    fq = Fourier_Quad(12, 124)
    bin_num2 = int(bin_num/2)
    inverse = range(int(bin_num / 2 - 1), -1, -1)
    mg_bin = fq.set_bin(mg, bin_num, 100)

    chisq_num = gh.shape[0]
    chisq = numpy.zeros((chisq_num,))
    for i in range(chisq_num):
        chisq[i] = fq.get_chisq(mg, mnu, gh[i], mg_bin, bin_num2, inverse, 0)
    return chisq


def img_text(ax, x, y, paras=None, strs_content=None, c="black", size=20):
    if paras:
        strs = "%.4f, %.5f, %.5f"%(paras[0], paras[1], paras[2])
        ax.text(x, y, strs, color=c, ha='left', va='center', transform=ax.transAxes, fontsize=size)
    if strs_content:
        ax.text(x, y, strs_content, color=c, ha='left', va='center', transform=ax.transAxes, fontsize=size)
