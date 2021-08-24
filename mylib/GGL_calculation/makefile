HK_FQLIB_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
HK_MPILIB_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_iolib.cpp -lhdf5  -lcurl
HK_GGL_FUNCTIONS_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_science_cal_lib.cpp


check:
	@echo $(HK_MYWORK_DIR)

CODE_NAME=GGL_cal.cpp
FINAL_NAME=GGL_cal

compile:
	mpicxx -I.  $(CODE_NAME)  $(HK_FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) $(HK_GGL_FUNCTIONS_PATH) -o $(FINAL_NAME) -std=c++11
clean:
	rm $(FINAL_NAME)


# CODE_NAME=test.cpp
# FINAL_NAME=test

# compile:
# 	g++ -I.  $(CODE_NAME)  $(FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) $(GGL_FUNCTIONS_PATH) -o $(FINAL_NAME) -std=c++11
# clean:
# 	rm $(FINAL_NAME)
