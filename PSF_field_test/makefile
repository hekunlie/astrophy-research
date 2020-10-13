FQLIB_PATH=/home/hklee/work/mylib/FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
HK_MPILIB_PATH=/home/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/home/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio -lcurl
COR_FUNCTIONS_APTH=/home/hklee/work/mylib/functions.cpp



CODE_NAME_S=simu.cpp
FINAL_NAME_S=simu

compile_s:
	mpicxx -I.  $(CODE_NAME_S) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) -o $(FINAL_NAME_S) -std=c++11
clean_s:
	rm $(FINAL_NAME_S)


CODE_NAME_C=calculate.cpp
FINAL_NAME_C=calculate

compile_c:
	mpicxx -I.  $(CODE_NAME_C) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) -o $(FINAL_NAME_C) -std=c++11
clean_c:
	rm $(FINAL_NAME_C)



