FQLIB_PATH=/home/hklee/work/mylib/FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
HK_MPILIB_PATH=/home/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/home/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio -lcurl
GGL_LIB_PATH=/home/hklee/work/mylib/GG_lensing_functions_expo_wise.cpp

# FQLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/FQlib.cpp -lgsl -lgslcblas -lfftw3l -lfftw3f -lm
# HK_MPILIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_mpi.cpp
# HK_IOLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio
# COR_FUNCTIONS_APTH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/functions.cpp
#GGL_LIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/GG_lensing_functions_expo_wise.cpp

CODE_NAME=GGL_cal.cpp
FINAL_NAME=GGL_cal

compile:
	g++ -I.  $(CODE_NAME) $(GGL_LIB_PATH)  -o $(FINAL_NAME) -std=c++11
clean:
	rm $(FINAL_NAME)


