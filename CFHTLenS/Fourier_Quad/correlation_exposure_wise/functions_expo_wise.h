#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include<hk_iolib.h>
#include<FQlib.h>
#include<vector>

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
    int mg1_idx;
    int mg2_idx;
    int mn_idx;
    int mu_idx;
    int mv_idx;
    int ra_idx;
    int dec_idx;
    int cos_dec_idx;
    int redshift_idx;

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
    // the original CFHT exposure label
    // to avoid the pairs in the same exposures
    int *obs_expo_label_1;
    int *obs_expo_label_2;

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

    int theta_accum_len_true;
    int expo_chi_block_len_true;

    // element num in buffer = 1024*1024*1024 Byte / 8 Byte, 8 Bytes for double, 
    // be carefull in some platforms of which the INT is short than 4 Bytes
    double *men_buffer;
    int max_buffer_size;
    int actual_buffer_size;
    int max_block_in_buffer;
    int block_size_in_buffer;
    int block_count;
    int *buffer_label;
    int total_buffer_num;
    int *task_expo_pair_jack_label;
    int jack_label_1, jack_label_2;



    double *corr_cal_stack_expo_theta_accum, *corr_cal_stack_expo_theta_num_accum;
    double *corr_cal_stack_num_count_chit, *corr_cal_stack_num_count_chix;

    // mpi task distribution
    int *expo_pair_num_each_rank;
    int my_expo_pair_st, my_expo_pair_ed;
    int *task_expo_label;
    int task_expo_num;
    std::vector<int>thread_pool;
    std::vector<int>thread_del;

    std::vector<int>task_expo_pair_labels_1;
    std::vector<int>task_expo_pair_labels_2;


    int result_file_tag;
    int task_complete;
    std::vector<int> expo_pair_label_1;
    std::vector<int> expo_pair_label_2;

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
};


void initialize(data_info *expo_info);

void line_count(char *file_path, data_info* expo_info);
void line_count(char *file_path, int &lines);

void read_list(char *file_path, data_info *expo_info, int &read_file_num);

void read_data(data_info *expo_info);

void read_expo_data_1(data_info *expo_info, int expo_label);
void read_expo_data_2(data_info *expo_info, int expo_label);

void initialize_expo_chi_block(data_info *expo_info);

void collect_chi_block(data_info *expo_info, int field_label);

void save_expo_data(data_info*field_info, int expo_label_1,int expo_label_2, int rank);
void save_expo_data_new(data_info *expo_info, int rank, int task_end_tag);

void save_expo_pair_label(data_info *expo_info, int rank);

void save_expo_data(data_info *expo_info, int expo_label, char *file_name);

void merge_data(data_info *expo_info);

void task_distribution(int portion, int my_id, data_info *field_info);

void task_prepare(int numprocs, int rank, data_info *field_info);

void initialize_thread_pool(data_info*expo_info,int numprocs);
void thread_pool_resize(data_info *expo_info);

void hist_2d_fast(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int bin_num, int bin_num2, int &ix, int &iy);
void hist_2d_new(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int bin_num, int bin_num1,int bin_num2, int bin_num3,int &ix, int &iy);
void hist_2d_new(MY_FLOAT x, MY_FLOAT y, MY_FLOAT*bins, int *bin_num_para,int &ix, int &iy);//faster than above

void hist_2d_new(MY_FLOAT*bins, int bin_num, MY_FLOAT *xy, int *bin_para, int &ix, int &iy);

void expo_distance(data_info *expo_info, int expo_label_0, int expo_label_1, int &label);
// if lable == 1, calculate, else, not

void find_pairs(data_info *field_info, int expo_label_0, int expo_label_1);
// read all exposures
void find_pairs_new(data_info *field_info, int expo_label_0, int expo_label_1, MY_FLOAT *gg1, MY_FLOAT *gg2);
// read the exposure needed
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



///////////////////////////////////// for the last step, \Xi^2 calculation and estimation of correlation function /////////////////////////////

struct corr_cal
{
    // for the last step of chi squared calculation, get_corr.cpp
    char parent_path[500];
    char result_path[550];
    char set_name[50];
    char log_path[550];
    char inform[600];
    
    // the guess of chi_{\pm} of PDF_SYM
    MY_FLOAT *chi_guess;
    int chi_guess_num;
    int chi_block_len, ir_chi_block_len, iz_chi_block_len,expo_chi_block_len;
    
    MY_FLOAT *mg_bin;
    int mg_bin_num, mg_bin_num1, mg_bin_num2, mg_bin_num3;

    int zbin_num;
    MY_FLOAT *zbin;
    /////////////////////////////////////

    // radius bin
    int theta_bin_num;
    MY_FLOAT *theta_bin;  

    double *expo_num_count_chit;
    double *expo_num_count_chix;
    double *theta_accum;
    double *theta_num_accum;

    int theta_accum_len_true;
    int expo_chi_block_len_true;
    int expo_block_len_in_buffer;
    
    double *corr_cal_expo_theta_accum[MAX_EXPO];
    double *corr_cal_expo_theta_num_accum[MAX_EXPO];
    double *corr_cal_expo_num_count_chit[MAX_EXPO];
    double *corr_cal_expo_num_count_chix[MAX_EXPO];
    
    double *corr_cal_stack_expo_theta_accum, *corr_cal_stack_expo_theta_num_accum;
    double *corr_cal_stack_num_count_chit, *corr_cal_stack_num_count_chix;
    
    int corr_cal_result_file_num;
    int *corr_cal_expo_pair_label[4];
    int *corr_cal_expo_pair_file_label;
    int corr_cal_total_pair_num;

    int corr_cal_expo_num;
    
    int expo_jack_id;
    int resample_num;
    int corr_cal_thread_num;
    int corr_cal_rank;
    int my_resample_label_st, my_resample_label_ed;

    // record the start & end expo label of each sub-sample
    // "resample_num" sub-samples
    int *jackknife_subsample_pair_st;
    int *jackknife_subsample_pair_ed;
    // for the MPI, distribute the resample tasks to "corr_cal_thread_num"
    // CPUs, these two arrays record the start&end resample-task label of each CPU
    int *jackknife_resample_st;
    int *jackknife_resample_ed;
    int jackknife_label;

    int corr_cal_chi_num;
    int corr_cal_final_data_num;
    double *corr_cal_mean_theta[MAX_RESAMPLE];
    double *corr_cal_chi_tt[MAX_RESAMPLE], *corr_cal_chi_xx[MAX_RESAMPLE];
    double *corr_cal_gtt[MAX_RESAMPLE], *corr_cal_gxx[MAX_RESAMPLE];
    double *corr_cal_gtt_sig[MAX_RESAMPLE], *corr_cal_gxx_sig[MAX_RESAMPLE];
    double *corr_cal_chi_guess;

};

void read_para(corr_cal *all_paras);

void prepare_data(corr_cal *all_paras, int tag);

void pre_jackknife(corr_cal *all_paras);

void resample_jackknife(corr_cal *all_paras,int resample_label);

void chisq_2d(double *num_count, int mg_bin_num, double &chisq);

void corr_calculate(corr_cal *all_paras, int resample_label);

void corr_task_alloc(int total_task_num, int portion, int *label_st, int *label_ed);

void save_result(corr_cal *all_paras);

#endif


