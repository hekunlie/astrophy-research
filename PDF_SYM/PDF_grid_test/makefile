# FQLIB_PATH=/home/hklee/work/mylib/FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
# HK_MPILIB_PATH=/home/hklee/work/mylib/hk_mpi.cpp
# HK_IOLIB_PATH=/home/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcurl


FQLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
HK_MPILIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcurl


CODE_NAME=source_simu.cpp
FINAL_NAME=simu

compile:
	mpicxx -I.  $(CODE_NAME) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) -o $(FINAL_NAME) -std=c++11
clean:
	rm $(FINAL_NAME)



# CODE_NAME=test.cpp
# FINAL_NAME=test

# compile_t:
# 	mpicxx -I.  $(CODE_NAME) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) -o $(FINAL_NAME) -std=c++11
# clean_t:
# 	rm $(FINAL_NAME)


