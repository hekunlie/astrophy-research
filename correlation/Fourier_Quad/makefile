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


CODE_NAME=calculate_correlation_expo_wise.cpp
FINAL_NAME=cal_cor

compile_expo:
	mpicxx -I.  $(CODE_NAME) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions_expo_wise.cpp -o $(FINAL_NAME) -std=c++11
clean_expo:
	rm $(FINAL_NAME)


CODE_NAME=calculate_correlation_expo_wise.cpp
FINAL_NAME=cal_cor_test

compile_t:
	mpicxx -I.  $(CODE_NAME) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions_expo_wise.cpp -o $(FINAL_NAME) -std=c++11
clean_t:
	rm $(FINAL_NAME)

CODE_NAME_t=speed_test.cpp
FINAL_NAME_t=test

compile_st:
	mpicxx -I.  $(CODE_NAME_t) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions.cpp -o $(FINAL_NAME_t) -std=c++11
clean_st:
	rm $(FINAL_NAME_t)

