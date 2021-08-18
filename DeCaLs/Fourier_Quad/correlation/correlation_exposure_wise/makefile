HK_FQLIB_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_FQlib.cpp  -lfftw3l -lfftw3f -lgsl -lgslcblas -lm
HK_MPILIB_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_mpi.cpp
HK_IOLIB_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_iolib.cpp -lhdf5  -lcurl
HK_COR_FUNCTIONS_PATH=$(HK_MYWORK_DIR)/work/mylib/hk_correlation_functions_expo_wise.cpp


check:
	@echo $(HK_MYWORK_DIR)

CODE_NAME_D=calculate_correlation_expo_wise_diff_expo.cpp
FINAL_NAME_D=cal_cor_diff_expo

CODE_NAME_S=calculate_correlation_expo_wise_same_expo.cpp
FINAL_NAME_S=cal_cor_same_expo

CODE_NAME_ST=calculate_correlation_expo_wise_stack_expo.cpp
FINAL_NAME_ST=cal_cor_stack_expo

compile_e:
	mpicxx -I.  $(CODE_NAME_D) $(HK_FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH)  $(HK_COR_FUNCTIONS_PATH) -O1 -o $(FINAL_NAME_D) -std=c++11
	mpicxx -I.  $(CODE_NAME_S) $(HK_FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH)  $(HK_COR_FUNCTIONS_PATH) -O1 -o $(FINAL_NAME_S) -std=c++11
	mpicxx -I.  $(CODE_NAME_ST) $(HK_FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH)  $(HK_COR_FUNCTIONS_PATH) -O1 -o $(FINAL_NAME_ST) -std=c++11

clean_e:
	rm $(FINAL_NAME_D)
	rm $(FINAL_NAME_S)
	rm $(FINAL_NAME_ST)



CODE_NAME_GC=get_corr.cpp
FINAL_NAME_GC=get_corr

compile_gc:
	mpicxx -I.  $(CODE_NAME_GC) $(HK_FQLIB_PATH) $(HK_MPILIB_PATH) $(HK_IOLIB_PATH) $(HK_COR_FUNCTIONS_PATH) -o $(FINAL_NAME_GC) -std=c++11
clean_gc:
	rm $(FINAL_NAME_GC)

