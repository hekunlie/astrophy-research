#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<hk_iolib.h>

#define MAX_FIELD 2000

struct data_info
{
    char *field_name[MAX_FIELD];

    ///////// for 2d correlation /////////
    float *field_data[MAX_FIELD];
    /////////////////////////////////////

    ////////// for tomography ///////////
    float *field_data_z1[MAX_FIELD];// field data, G1, G2 .. of zbin 1
    float *field_data_z2[MAX_FIELD];// field data, G1, G2 .. of zbin 2
    int *total_gal_num_z1;// gal num in each field of zbin 1
    int *total_gal_num_z2;// gal num in each field of zbin 2
    // redshift bin
    int zbin_num;
    int *num_in_zbin;
    float *zbin;
    /////////////////////////////////////

    // radius bin
    int theta_bin_num;
    float *theta_bin;

    int field_data_col;

    int total_field_num;

    int *exposure_num_of_field;
    
    float *field_cen_ra;
    float *field_cen_dec;
    float *field_cen_cos_dec;
    float *delta_ra;
    float *delta_dec;
    float *delta_len;
    
    // mpi task distribution
    int *field_num_each_rank;
    int my_field_st, my_field_ed;
};


void read_file(char *file_path, data_info *field_info, int &read_file_num);

void read_field_data(data_info *field_info, int zbin_label_0, int zbin_label_1);

void initialize(char *file_path, data_info *field_info, int total_field_num, int numprocs, int rank);

void task_distribution(int portion, int my_id, data_info *field_info);

void fast_hist(float data, float*bins, int *num_in_bin, int bin_num);

void find_pairs(data_info *field_info);
#endif


