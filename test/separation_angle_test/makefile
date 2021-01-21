FQLIB_PATH=/home/hklee/work/mylib/FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
HK_MPILIB_PATH=/home/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/home/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio -lcurl
SCI_LIB_PATH=/home/hklee/work/mylib/science_cal_lib.cpp

# FQLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/FQlib.cpp -lgsl -lgslcblas -lfftw3l -lfftw3f -lm
# HK_MPILIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_mpi.cpp
# HK_IOLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio
# SCI_LIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/GG_lensing_functions_expo_wise.cpp

CODE_NAME=separation_angle_test.cpp
FINAL_NAME=separation_angle_test

compile:
	mpicxx -I.  $(CODE_NAME)  $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) $(SCI_LIB_PATH) -o $(FINAL_NAME) -std=c++11
clean:
	rm $(FINAL_NAME)


# CODE_NAME=test.cpp
# FINAL_NAME=test

# compile:
# 	g++ -I.  $(CODE_NAME)  $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) $(SCI_LIB_PATH) -o $(FINAL_NAME) -std=c++11
# clean:
# 	rm $(FINAL_NAME)
