# FQLIB_PATH=/home/hklee/work/mylib/FQlib.cpp -lgsl -lgslcblas -lfftw3l -lfftw3f -lm
# HK_MPILIB_PATH=/home/hklee/work/mylib/hk_mpi.cpp
# HK_IOLIB_PATH=/home/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio
# COR_FUNCTIONS_APTH=/home/hklee/work/mylib/functions.cpp

FQLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/FQlib.cpp -lgsl -lgslcblas -lfftw3l -lfftw3f -lm
HK_MPILIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio
COR_FUNCTIONS_APTH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/functions.cpp


CODE_NAME_F=calculate_correlation.cpp
FINAL_NAME_F=cal_cor

compile_f:
	mpicxx -I.  $(CODE_NAME_F) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions.cpp -o $(FINAL_NAME_F) -std=c++11
clean_f:
	rm $(FINAL_NAME_F)


CODE_NAME_E=calculate_correlation_expo_wise.cpp
FINAL_NAME_E=cal_cor

compile_e:
	mpicxx -I.  $(CODE_NAME_E) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions_expo_wise.cpp -O -o $(FINAL_NAME_E) -std=c++11
clean_e:
	rm $(FINAL_NAME_E)


CODE_NAME_T=calculate_correlation_expo_wise.cpp
FINAL_NAME_T=cal_cor_test

compile_t:
	mpicxx -I.  $(CODE_NAME_T) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions_expo_wise.cpp -o $(FINAL_NAME_T) -std=c++11
clean_t:
	rm $(FINAL_NAME_T)


CODE_NAME_TO=calculate_correlation_expo_wise.cpp
FINAL_NAME_TO=cal_cor_test_O

compile_to:
	mpicxx -I.  $(CODE_NAME_TO) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions_expo_wise.cpp -O -o $(FINAL_NAME_TO) -std=c++11
clean_to:
	rm $(FINAL_NAME_TO)


CODE_NAME_st=speed_test.cpp
FINAL_NAME_st=test

compile_st:
	mpicxx -I.  $(CODE_NAME_st) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) functions.cpp -o $(FINAL_NAME_st) -std=c++11
clean_st:
	rm $(FINAL_NAME_st)

