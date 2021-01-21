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




const double C_0_hat = 2.99792458; // 10^8
const double C_0 = 2.99792458*1.e8;// speed of light m/s

const double H_0_hat = 0.70;
const double H_0 = 70; //  Km/s/Mpc

const double G_0_hat= 6.6740831; //  10^{-11}
const double G_0 = 6.6740831*1.e-11; //  m^3 s^{-2} Kg^{-1}

const double One_Light_Year_hat = 9.4607304725808;// 10^15
const double One_Light_Year = 9.4607304725808*1.e15;// meter

const double One_Mpc_hat = 3.085677581; // 10^22
const double One_Mpc = 3.085677581*1.e22;// meter

const double M_sun_hat = 1.9885;
const double M_sun = 1.9885*1.e30;//Kg

const double DEG2RAD = 0.01745329251994329576923;// pi/180

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
    int mg_bin_num;
    int pdf_guess_num;

    double *gt_guess;
    double *delta_sigma_guess;
    MY_FLOAT *mg_bin;

    int chi_theta_block_len, chi_jack_block_len;
    int chi_signal_block_len;
    int total_chi_count_len;

    // total source count in the PDF, 7 parts
    // the first one, len = chi_signal_block_len, for the \Delta\Sigma, excess surface density
    // the second one, len = chi_signal_block_len, for the cross \Delta\Sigma, should be consitent with 0, systematic check
    // the third one, len = chi_signal_block_len, for the tangential \gamma,
    // the forth one, len = chi_signal_block_len, for the cross \gamma, should be consitent with 0, systematic check

    // the fifth one, len = signal_pts_num, for theta, separation angle, radian
    // the sixth one, len = signal_pts_num, for comoving distance, separation comoving distance, Mpc
    // the seventh one, len = signal_pts_num, for pair count,
    double *total_chi_count;        
    
    // for each individual calculation, it will be added to the total one when finished
    double *worker_sub_chi_count;
    double *worker_total_chi_count;

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
    int len_prop_dist_col;
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
    int src_prop_dist_col;
    

    //////////////////  task distribution  //////////////////////
    std::vector<int> task_len_expo_labels;
    std::vector<int> task_src_expo_labels;
    int task_expo_num;
};

void ggl_initialize(ggl_data_info *data_info);

// void ggl_task_prepare(ggl_data_info *data_info);

void ggl_read_len_list(char *file_path, ggl_data_info* expo_info);
void ggl_read_src_list(char *file_path, ggl_data_info* expo_info);

void ggl_read_len_exp(ggl_data_info *data_info, int len_expo_label);
void ggl_read_src_exp(ggl_data_info *data_info, int src_expo_label);

void ggl_find_src_needed(ggl_data_info *data_info, int len_expo_label);

void ggl_rotation_matrix(MY_FLOAT cent_ra, MY_FLOAT cent_dec, MY_FLOAT cent_cos_dec, MY_FLOAT src_ra, MY_FLOAT src_dec, MY_FLOAT*rotation_matrix);

void ggl_rotate_estimator(MY_FLOAT G1, MY_FLOAT G2, MY_FLOAT U, MY_FLOAT V, MY_FLOAT *rotation_matrix, MY_FLOAT &Gt, MY_FLOAT &Gx, MY_FLOAT &Ut);

void ggl_fast_hist(MY_FLOAT *bins, int bin_num, MY_FLOAT val, int pre_bin_tag, int &bin_tag);

void ggl_find_pair(ggl_data_info *data_info, int len_expo_label);

void ggl_collect_chi(ggl_data_info *data_info);

void ggl_cal_signals(ggl_data_info * data_info);

void ggl_pdf_signals(double *chi_count, double*pdf_signal_guess, int pdf_guess_num, int mg_bin_num, int signal_pts_num, double *signal, double *signal_err);

#endif
