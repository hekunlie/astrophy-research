#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<hk_iolib.h>

#define MAX_FIELD 2000
#define MAX_EXPOS_NUM 20000


struct data_info
{
    char *field_name[MAX_FIELD];
    char *exposure_name[MAX_EXPOS_NUM];

    int total_exposure_num;

    int *field_label;
    int *exposure_label;
    int *exposure_num_of_field;
    
    float *field_cen_ra;
    float *field_cen_dec;
    float *delta_ra;
    float *delta_dec;

    // mpi task distribution
    int my_exposure_st;
    int my_exposure_ed;
    int my_exposure_num;
};


void read_file(char *file_path, data_info *field_info, int &read_file_num);

void initialize(char *file_path, data_info *field_info, int total_exposure_num, int numprocs, int rank);
// allocate the array & read the file of information of all exposures

void task_distribution(int total_task_num, int portion, int my_id, data_info *field_info);


void find_pairs(data_info *field_info);
#endif


