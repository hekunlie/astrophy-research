from sys import path
path.append("D:/Github/astrophy-research/mylib")
path.append("D:/Github/astrophy-research/multi_shear_detect")
from Fourier_Quad import Fourier_Quad
import numpy

# whatever number in the bracket,
# is not for shear measurement, just a small one is ok
fq = Fourier_Quad(1, 123)

# read data from catalog
# the standard output files are "***.cat" of each exposure
# please read the its header before any work

################# source selection ####################
# it provides selection criteria for source selection

# nstar, shoudl be >= 12 or larger, is the star number on the chip,
# the PSF reconstruction will fail if the chip has too less stars,
# in the latest version, we use >= 20

# flux2 or flux_alt, is our SNR defined in Fourier Space
# >= 3 or 4 is fine, you can test it

# gf1, gf2 are the field distortion, they are actually the components
# of ellipticity induced by coordinates projection
# <= 0.005 for both is safe

# if you have some other selection criteria, it is ok to use it simultaneously
# (*I think, you can test it*)

################# shear measurement ####################
# g1, g2, de, h1, h2 are needed for shear measurement
# they are just the G1, G2, N, U, and V defined in the paper and Fourier_Quad class

# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# !!!! remember put a '-' to h1 and h2 before calculation
# !!!! because the definition in Jun's pipeline is opposite to the one defined in the paper and my code
h1 = -h1
h2 = -h2
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!remember it

# the following part is the calculation
mg1 = g1
mg2 = g2
mn = de
mu = h1
mv = h2
# if you need rotate the galaxy or shear estimators like ellipticity rotation, you can do as:
# !!! 2theta
mg1 = numpy.cos(2*theta)*g1 - numpy.sin(2*theta)*g2
mg2 = numpy.sin(2*theta)*g1 + numpy.cos(2*theta)*g2
# N (de) is a scalar
mn = de
# !!! 4theta
mu = numpy.cos(4*theta)*h1 - numpy.sin(4*theta)*h2
mv = numpy.sin(4*theta)*h1 + numpy.cos(4*theta)*h2

# for g1 (the shear) estimating, call the PDF_SYM method in Fourier_Quad
g1, g1_sig = fq.find_shear(mg1, mn + mu, bin_num=8)[:2]
# for g2 (the shear) estimating
g2, g2_sig = fq.find_shear(mg2, mn - mu, bin_num=8)[:2]
# you can use large bin number for PDF_SYM, say 10 or 12. However, 8 is enough.
# If you do not have enough galaxies, say a few hundreds, you shoudl use 4 bins for PDF_SYM
