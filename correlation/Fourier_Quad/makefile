FQLIB_PATH=/home/hklee/work/mylib/FQlib.cpp -lgsl -lgslcblas -lfftw3l -lfftw3f -lm
HK_MPILIB_PATH=/home/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/home/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio
COR_FUNCTIONS_APTH=/home/hklee/work/mylib/functions.cpp

CODE_NAME=calculate_correlation.cpp
FINAL_NAME=cal_cor

compile:
	mpicxx -I.  $(CODE_NAME) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions.cpp -o $(FINAL_NAME) -std=c++11
clean:
	rm $(FINAL_NAME)
