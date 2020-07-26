FQLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/FQlib.cpp -lgsl -lgslcblas -lfftw3l -lfftw3f -lm
HK_MPILIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=/lustre/home/acct-phyzj/phyzj-sirius/hklee/work/mylib/hk_iolib.cpp -lhdf5 -lcfitsio

CODE_NAME_m=check_ct_property.cpp
FINAL_NAME_m=check_ct
compile:
	mpicxx -I.  $(CODE_NAME_m) $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) -o $(FINAL_NAME_m) -std=c++11
clean:
	rm $(FINAL_NAME_m)

