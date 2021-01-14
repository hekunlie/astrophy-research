#ifndef GGL_FUNCTIONS_H
#define GGL_FUNCTIONS_H

#include<hk_iolib.h>
#include<FQlib.h>


#define MAX_JACK 2000
#define MAX_EXPO_NUM 50000
#define MY_FLOAT float
#define GGL_PROP_DIST_STACK

struct ggl_data_info
{
    int jack_id;
    int jack_num;
    // the signal number, theta or comoving distance point number
    int signal_pts_num;
    
    MY_FLOAT *sepration_bin;

    int pair_count;

    ////////////////////  for the SYM_PDF method of Fourier_Quad  ///////////////////////
    // bin num for G1,G2
    int mg_bin_num;
    int pdf_guess_num;

    MY_FLOAT *gt_guess;
    MY_FLOAT *delta_sigma_guess;
    MY_FLOAT *mg_bin;

    int chi_len;
    int signal_chi_len;


    // total source count in the PDF, 7 parts
    // the first one, len = signal_chi_len, for the \Delta\Sigma, excess surface density
    // the second one, len = signal_chi_len, for the cross \Delta\Sigma, should be consitent with 0, systematic check
    // the third one, len = signal_chi_len, for the tangential \gamma,
    // the forth one, len = signal_chi_len, for the cross \gamma, should be consitent with 0, systematic check

    // the fifth one, len = signal_pts_num, for theta, separation angle, radian
    // the sixth one, len = signal_pts_num, for comoving distance, separation comoving distance, Mpc
    // the seventh one, len = signal_pts_num, for pair count,
    double *total_chi_count;        
    
    // for each individual calculation, it will be added to the total one when finished
    double *indi_chi_count;


    int pos_inform_num;

    ///////////////////  the informs of each len exposure //////////////
    int len_expo_label;
    // the position informs of each len exposure file
    char *len_expo_path[MAX_EXPO_NUM];
    int *len_data_row;
    int len_data_col;
    int len_expo_num;
  
    MY_FLOAT *len_expo_data;
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


    MY_FLOAT back_dz;
    MY_FLOAT *separation_bin;
    int sep_bin_num;
    
};

void ggl_initialize(ggl_data_info *data_info);

void ggl_read_len_exp(ggl_data_info *data_info, int len_expo_label);
void ggl_read_src_exp(ggl_data_info *data_info, int src_expo_label);

void ggl_find_src_needed(ggl_data_info *data_info, int len_expo_label);

#endif
