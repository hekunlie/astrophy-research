#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<hk_iolib.h>
#include<FQlib.h>

#define MAX_EXPO 2000
#define MAX_RESAMPLE 2000
#define MY_FLOAT float

struct data_info
{
    char parent_path[500];
    char *expo_name_path[MAX_EXPO];
    char *expo_name[MAX_EXPO];
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

    int total_expo_num;

    MY_FLOAT *expo_data[MAX_EXPO];
    int expo_data_col;

    int *expo_num_in_zbin[MAX_EXPO];// gal num in each zbin of each exposure
    int *expo_zbin_st[MAX_EXPO];// the start & end row of each zbin in each exposure
    int *expo_zbin_ed[MAX_EXPO];
    int *expo_zbin_label[MAX_EXPO];// label of zbin, each gal

    int data_read_label_1, data_read_label_2;
    MY_FLOAT *expo_data_1;
    MY_FLOAT *expo_data_2;
    int *expo_zbin_label_1;
    int *expo_zbin_label_2;

    int *expo_gal_num;// gal num in each exposure
    // about the position of each exposure
    MY_FLOAT *expo_cen_ra;
    MY_FLOAT *expo_cen_dec;
    MY_FLOAT *expo_cen_cos_dec;
    MY_FLOAT *expo_delta_ra; // half width of RA
    MY_FLOAT *expo_delta_dec;// half width of DEC
    MY_FLOAT *expo_delta_len; // sqrt((delta_ra*cen_cos_dec)^2 + delta_dec^2)
    // redshift bin
    int zbin_num;
    MY_FLOAT *zbin;
    /////////////////////////////////////

    // radius bin
    int theta_bin_num;
    MY_FLOAT *theta_bin;  
    // save the separation of each pair
    MY_FLOAT *theta;
    // for the calculation of the mean theta at last
    // the x position in the last figure
    double *theta_accum, *theta_num_accum;
    int theta_accum_len;


    // mpi task distribution
    int *expo_pair_num_each_rank;
    int my_expo_pair_st, my_expo_pair_ed;
    int *task_expo_label;
    int task_expo_num;


    // the guess of chi_{\pm} of PDF_SYM
    MY_FLOAT *chi_guess;
    int chi_guess_num;
    int chi_block_len, ir_chi_block_len, iz_chi_block_len,expo_chi_block_len;
    
    MY_FLOAT *mg_bin;
    int mg_bin_num, mg_bin_num1, mg_bin_num2, mg_bin_num3;


    // of each field, for errorbar estimation
    double *num_count_chit[MAX_EXPO];
    double *num_count_chix[MAX_EXPO];
    // number counting of each exposure for the signal estimation
    double *expo_num_count_chit;
    double *expo_num_count_chix;
    double gg_pairs;
    
    int gg_len;
    MY_FLOAT *gg_1;
    MY_FLOAT *gg_2;
    int loop_label;

    // for the last step of chi squared calculation, get_corr.cpp
    double *corr_cal_expo_theta_accum[MAX_EXPO];
    double *corr_cal_expo_theta_num_accum[MAX_EXPO];
    double *corr_cal_expo_num_count_chit[MAX_EXPO];
    double *corr_cal_expo_num_count_chix[MAX_EXPO];
    
    double *corr_cal_stack_expo_theta_accum, *corr_cal_stack_expo_theta_num_accum;
    double *corr_cal_stack_num_count_chit, *corr_cal_stack_num_count_chix;
    
    int theta_accum_len_true;
    int expo_chi_block_len_true;

    int resample_num;
    int corr_cal_chi_len;
    int corr_cal_final_data_num;
    double *corr_cal_mean_theta;
    double *corr_cal_chi_tt[MAX_RESAMPLE], *corr_cal_chi_xx[MAX_RESAMPLE];
    double *corr_cal_gtt[MAX_RESAMPLE], *corr_cal_gxx[MAX_RESAMPLE];
    double *corr_cal_gtt_sig[MAX_RESAMPLE], *corr_cal_gxx_sig[MAX_RESAMPLE];
};


void initialize(data_info *field_info, int total_field_num);

void read_list(char *file_path, data_info *field_info, int &read_file_num);

void read_data(data_info *field_info);

void read_expo_data_1(data_info *field_info, int expo_label);
void read_expo_data_2(data_info *field_info, int expo_label);

void initialize_expo_chi_block(data_info *field_info);

void collect_chi_block(data_info *field_info, int field_label);

void save_expo_data(data_info*field_info, int expo_label);

void save_expo_data(data_info *expo_info, int expo_label, char *file_name);

void task_distribution(int portion, int my_id, data_info *field_info);

void task_prepare(int numprocs, int rank, data_info *field_info);

void hist_2d_fast(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int bin_num, int bin_num2, int &ix, int &iy);
void hist_2d_new(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int bin_num, int bin_num1,int bin_num2, int bin_num3,int &ix, int &iy);
void hist_2d_new(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int *bin_num_para,int &ix, int &iy);//faster than above

void hist_2d_new(MY_FLOAT*bins, int bin_num, MY_FLOAT *xy, int *bin_para, int &ix, int &iy);

void expo_distance(data_info *expo_info, int expo_label_0, int expo_label_1, int &label);
// if lable == 1, calculate, else, not

void find_pairs(data_info *field_info, int expo_label_0, int expo_label_1);
// read all exposures
void find_pairs_new(data_info *field_info, int expo_label_0, int expo_label_1);
// read the exposure needed
#endif


