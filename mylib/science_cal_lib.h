#ifndef SCI_LIB_H
#define SCI_LIB_H


#include<hk_iolib.h>
#include<hk_mpi.h>
#include<FQlib.h>
#include<vector>


#define MAX_JACK 2000
#define MAX_EXPO_NUM 50000
#define MY_FLOAT float
#define MY_MPI_TYPE MPI_FLOAT
// #define GGL_PROP_DIST_STACK
#define GGL_COM_DIST_STACK
#define GGL_GAMMA_T
#define GGL_DELTA_SIGMA

const double C_0_hat = 2.99792458; //
const double C_0_scale = 1.e8;// speed of light m/s

const double G_0_hat= 6.6740831; // 
const double G_0_scale = 1.e-11; //  m^3 s^{-2} Kg^{-1}

const double One_Light_Year_hat = 9.4607304725808;// 10^15
const double One_Light_Year_scale = 1.e15;// meter

const double One_pc_hat = 3.085677581;
const double One_pc_scale = 1.e16;//meter

const double M_sun_hat = 1.9885;
const double M_sun_scale = 1.e30;//Kg

const double SCI_PI = 3.141592653589793;
const double DEG2RAD = 0.017453292519943;// pi/180
const double RAD2DEG = 57.295779513082323;


void separation_angle_1(const double RA1, const double DEC1, const double RA2, const double DEC2, double &sep_radian);
void separation_angle_1(const float RA1, const float DEC1, const float RA2, const float DEC2, float &sep_radian);

void separation_angle_2(const double RA1, const double DEC1, const double RA2, const double DEC2, double &sep_radian);// checked
void separation_angle_2(const float RA1, const float DEC1, const float RA2, const float DEC2, float &sep_radian);
/* calculate the separation angle between two points on the sphere */
/* RA & DEC in unit of degree */
/* https://en.wikipedia.org/wiki/Great-circle_distance */


void com_distance(const double low_z, const double high_z, const double omg_m, const double omg_lam, double &result, const double precision_thresh, const bool integ_only);


/////////////////////////  GGL part  /////////////////////////////////////////
struct ggl_data_info
{   
    char ggl_total_path[400];
    char set_name[50];
    char ggl_pdf_inform_path[450];
    char ggl_foreground_inform_path[450];
    char ggl_background_inform_path[450];

    char ggl_log_inform[500];
    char ggl_log_path[500];

    char ggl_result_path[500];

    int jack_id;
    int jack_num;
    
    // the signal number, theta or comoving distance point number
    int signal_pts_num;
    
    MY_FLOAT back_dz;
    MY_FLOAT *separation_bin;
    int sep_bin_num;

    int pair_count;

    int rank, numprocs;


    ////////////////////  for the SYM_PDF method of Fourier_Quad  ///////////////////////
    // bin num for G1,G2
    int mg_sigma_bin_num;
    int mg_gt_bin_num;
    MY_FLOAT *mg_sigma_bin;
    MY_FLOAT *mg_gt_bin;

    int pdf_gt_num, pdf_sigma_num;
    double *gt_guess;
    double *delta_sigma_guess;

    int chi_sigma_theta_block_len, chi_sigma_theta_block_len_sub;
    int chi_sigma_jack_block_len;
    int chi_g_theta_block_len, chi_g_theta_block_len_sub;
    int chi_g_jack_block_len;

    int sub_signal_count_len, total_signal_count_len;

    double *total_chi_sigma_tan;
    double *total_chi_sigma_cross;
    double *total_chi_g_tan;        
    double *total_chi_g_cross;        
    double *total_signal_count;
    
    // for each individual calculation, it will be added to the total one when finished
    double *worker_sub_chi_sigma_tan;
    double *worker_total_chi_sigma_tan;

    double *worker_sub_chi_sigma_cross;
    double *worker_total_chi_sigma_cross;

    double *worker_sub_chi_g_tan;
    double *worker_total_chi_g_tan;
    double *worker_sub_chi_g_cross;
    double *worker_total_chi_g_cross;

    double *worker_sub_signal_count;
    double *worker_total_signal_count;


    MY_FLOAT crit_coeff;
    
    int pos_inform_num;

    ///////////////////  the informs of each len exposure //////////////
    int len_expo_label;
    // the position informs of each len exposure file
    
    char *len_expo_path[MAX_EXPO_NUM];
    char *len_expo_name[MAX_EXPO_NUM];
    int *len_data_row;
    int *len_expo_jackid;
    int len_data_col;
    int len_expo_num;
  
    MY_FLOAT *len_expo_data;
    MY_FLOAT *len_width_informs; // the width of each len exposure file in unit of radian
    MY_FLOAT *len_nearest_dist; // comoving or physical distance
    int len_expo_read_tag;

    int len_ra_col;
    int len_dec_col;
    int len_cos_dec_col;
    int len_z_col;
    int len_com_dist_col;
    int len_jackid_col;


    ///////////////////  the informs of each source exposure //////////////
    MY_FLOAT *src_pos_informs[MAX_EXPO_NUM];
    char *src_expo_path[MAX_EXPO_NUM];
    char *src_expo_name[MAX_EXPO_NUM];
    int *src_data_row;
    int *src_expo_needed_tag;
    int src_data_col;
    int src_expo_num;

    MY_FLOAT *src_expo_data;
    int src_expo_read_tag;
    int src_mg1_col;
    int src_mg2_col;
    int src_mn_col;
    int src_mu_col;
    int src_mv_col;

    int src_ra_col;
    int src_dec_col;
    int src_cos_dec_col;
    int src_z_col;
    int src_zerr_col;
    int src_com_dist_col;
    

    //////////////////  task distribution  //////////////////////
    std::vector<int> task_len_expo_labels;
    std::vector<int> task_src_expo_labels;
    int task_expo_num;
};

void ggl_initialize(ggl_data_info *data_info);

// void ggl_task_prepare(ggl_data_info *data_info);

void ggl_read_list(ggl_data_info *data_info);
void ggl_read_len_list(char *file_path, ggl_data_info* data_info);
void ggl_read_src_list(char *file_path, ggl_data_info* data_info);

void ggl_read_pdf_inform(ggl_data_info *data_info);

void ggl_read_len_exp(ggl_data_info *data_info, int len_expo_label);
void ggl_read_src_exp(ggl_data_info *data_info, int src_expo_label);

void ggl_set_bins(ggl_data_info *data_info);

void ggl_find_src_needed(ggl_data_info *data_info, int len_expo_label);

void ggl_rotation_matrix(MY_FLOAT cent_ra, MY_FLOAT cent_dec, MY_FLOAT cent_cos_dec, MY_FLOAT src_ra, MY_FLOAT src_dec, MY_FLOAT*rotation_matrix);

void ggl_rotate_estimator(MY_FLOAT G1, MY_FLOAT G2, MY_FLOAT U, MY_FLOAT V, MY_FLOAT *rotation_matrix, MY_FLOAT &Gt, MY_FLOAT &Gx, MY_FLOAT &Ut);

void ggl_fast_hist(MY_FLOAT *bins, int bin_num, MY_FLOAT val, int pre_bin_tag, int &bin_tag);

void ggl_find_pair(ggl_data_info *data_info, int len_expo_label);

void ggl_collect_chi(ggl_data_info *data_info);

void ggl_cache(ggl_data_info *data_info);

void ggl_cal_signals(ggl_data_info * data_info);

void ggl_pdf_signals(double *chi_count, double*pdf_signal_guess, int pdf_guess_num, int mg_bin_num, int signal_pts_num, double *signal, double *signal_err);

#endif
