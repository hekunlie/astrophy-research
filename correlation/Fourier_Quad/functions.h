#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<hk_iolib.h>

#define MAX_FIELD 2000

struct data_info
{
    char *field_name[MAX_FIELD];
    
    int total_field_num;

    int *exposure_num_of_field;
    
    float *field_cen_ra;
    float *field_cen_dec;
    float *delta_ra;
    float *delta_dec;
    float *delta_len;
    
    // mpi task distribution
    int *field_num_each_rank;
    int my_field_st, my_field_ed;
};


void read_file(char *file_path, data_info *field_info, int &read_file_num);

void initialize(char *file_path, data_info *field_info, int total_field_num, int numprocs, int rank);
// allocate the array & read the file of information of all exposures

void task_distribution(int portion, int my_id, data_info *field_info);


void find_pairs(data_info *field_info);
#endif


