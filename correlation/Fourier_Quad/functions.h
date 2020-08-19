#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<hk_iolib.h>
#include<FQlib.h>

#define MAX_FIELD 2000

#define MY_FLOAT float

struct data_info
{
    char parent_path[500];
    char *field_name_path[MAX_FIELD];
    char *field_name[MAX_FIELD];
    // data column index meaning (defined in the prepare_data.py)
    // remind to check the index before running
    int mg1_idx = 0;
    int mg2_idx = 1;
    int mn_idx = 2;
    int mu_idx = 3;
    int mv_idx = 4;
    int ra_idx = 5;
    int dec_idx = 6;
    int cos_dec_idx = 7;
    int redshift_idx = 8;

    ///////// for 2d correlation /////////
    MY_FLOAT *field_data[MAX_FIELD];
    /////////////////////////////////////

    ////////// for tomography ///////////
    int zbin_label_0, zbin_label_1;
    MY_FLOAT *field_data_z1[MAX_FIELD];// field data, G1, G2 .. of zbin 1
    MY_FLOAT *field_data_z2[MAX_FIELD];// field data, G1, G2 .. of zbin 2
    
    int *field_block_label_z1[MAX_FIELD]; // block labels
    int *field_block_label_z2[MAX_FIELD];
    int *field_expo_label_z1[MAX_FIELD]; // exposure labels
    int *field_expo_label_z2[MAX_FIELD];

    int *total_gal_num_z1;// gal num in each field of zbin 1
    int *total_gal_num_z2;// gal num in each field of zbin 2
    // redshift bin
    int zbin_num;
    int *num_in_zbin;
    MY_FLOAT *zbin;
    /////////////////////////////////////

    // radius bin
    int theta_bin_num;
    MY_FLOAT *theta_bin;

    int field_data_col;

    int total_field_num;

    int *exposure_num_of_field;
    
    // about the field
    MY_FLOAT *field_cen_ra;
    MY_FLOAT *field_cen_dec;
    MY_FLOAT *field_cen_cos_dec;
    MY_FLOAT *field_delta_ra;
    MY_FLOAT *field_delta_dec;
    MY_FLOAT *field_delta_len;
    
    // about the block in the field
    MY_FLOAT *block_st_z1[MAX_FIELD];
    MY_FLOAT *block_st_z2[MAX_FIELD];
    MY_FLOAT *block_ed_z1[MAX_FIELD];
    MY_FLOAT *block_ed_z2[MAX_FIELD];

    MY_FLOAT *block_cen_ra[MAX_FIELD];
    MY_FLOAT *block_cen_dec[MAX_FIELD];
    MY_FLOAT *block_cen_cos_dec[MAX_FIELD];
    MY_FLOAT *block_delta_len[MAX_FIELD];
    MY_FLOAT block_size;
    int *block_num;

    // mpi task distribution
    int *field_num_each_rank;
    int my_field_st, my_field_ed;

    // the guess of chi_{\pm} of PDF_SYM
    MY_FLOAT *chi_guess;
    int chi_guess_num;

    int chi_bin_num=10;
    double *num_count_chit[MAX_FIELD];
    double *num_count_chix[MAX_FIELD];

    int gg_len;
    MY_FLOAT *gg_tt[200];
    MY_FLOAT *gg_xx[200];

};


void read_inform(char *file_path, data_info *field_info, int &read_file_num);

void read_field_data(data_info *field_info);

void initialize(char *file_path, data_info *field_info, int total_field_num, int numprocs, int rank);

void task_distribution(int portion, int my_id, data_info *field_info);

void fast_hist(MY_FLOAT data, MY_FLOAT*bins, int *num_in_bin, int bin_num);

void find_pairs_same_field(data_info *field_info, int field_label);
void find_pairs_diff_field(data_info *field_info, int field_label_0, int field_label_1);

#endif


