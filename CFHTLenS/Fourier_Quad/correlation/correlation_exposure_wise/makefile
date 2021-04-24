FQLIB_PATH=/home/hklee/work/mylib/FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
HK_MPILIB_PATH=/home/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/home/hklee/work/mylib/hk_iolib.cpp -lhdf5  -l curl
COR_FUNCTIONS_PATH=/home/hklee/work/mylib/correlation_functions_expo_wise.cpp

# FQLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/FQlib.cpp -lgsl -lgslcblas -lfftw3l -lfftw3f -lm
# HK_MPILIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_mpi.cpp
# HK_IOLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_iolib.cpp -lhdf5 
# COR_FUNCTIONS_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/correlation_functions_expo_wise.cpp


CODE_NAME_E=calculate_correlation_expo_wise.cpp
FINAL_NAME_E=cal_cor

compile_e:
	mpicxx -I.  $(CODE_NAME_E) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH)  $(COR_FUNCTIONS_PATH) -O1 -o $(FINAL_NAME_E) -std=c++11
clean_e:
	rm $(FINAL_NAME_E)



CODE_NAME_GC=get_corr.cpp
FINAL_NAME_GC=get_corr

compile_gc:
	mpicxx -I.  $(CODE_NAME_GC) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) $(COR_FUNCTIONS_PATH) -o $(FINAL_NAME_GC) -std=c++11
clean_gc:
	rm $(FINAL_NAME_GC)

