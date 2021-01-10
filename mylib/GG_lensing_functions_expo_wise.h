#ifndef GGLEN_FUNCTIONS_H
#define GGLEN_FUNCTIONS_H

#define MAX_JACK 2000

#define MY_FLOAT float

struct ggl_data_info
{
    int jack_id;
    int jack_num;
    // the signal number, theta or comoving distance point number
    int signal_pts_num;
    
    int pair_count;
    int foreground_expo_label;
    int background_expo_label;

    ////////////////////  for the SYM_PDF method of Fourier_Quad  ///////////////////////
    // bin num for G1,G2
    int mg_bin_num;
    int pdf_guess_num;

    int chi_len;
    int signal_chi_len;
    // equals the signal number, theta or comoving distance point number
    int signal_pts_chi_len;
    
    // total source count in the PDF, 7 parts
    // the first one, len = signal_chi_len, for the \Delta\Sigma, excess surface density
    // the second one, len = signal_chi_len, for the cross \Delta\Sigma, should be consitent with 0, systematic check
    // the third one, len = signal_chi_len, for the tangential \gamma,
    // the forth one, len = signal_chi_len, for the cross \gamma, should be consitent with 0, systematic check

    // the fifth one, len = signal_pts_chi_len, for theta, separation angle, radian
    // the sixth one, len = signal_pts_chi_len, for comoving distance, separation comoving distance, Mpc
    // the seventh one, len = signal_pts_chi_len, for pair count,
    double *total_chi_count;        
    
    // for each individual calculation, it will be added to the total one when finished
    double *indi_chi_count;


    int pos_inform_num;

    ///////////////////  the informs of each foreground exposure //////////////
    // the position informs of each foreground exposure file
    MY_FLOAT *all_foreground_pos_infoms;
    int foreground_expo_num;
    int pos_inform_num;

    MY_FLOAT *foreground_expo_data;
    int foreground_data_col;
    int foreground_data_row;
 
    int foreground_ra_col;
    int foreground_dec_col;
    int foreground_cos_dec_col;
    int foreground_z_col;
    int foreground_com_dist_col;


    ///////////////////  the informs of each background exposure //////////////
    MY_FLOAT *all_background_informs;
    int background_data_num;
    int background_data_col;

    int background_ra_col;
    int background_dec_col;
    int background_cos_dec_col;
    int background_z_col;
    int background_com_dist_col;

    int background_mg1_col;
    int background_mg2_col;
    int background_mn_col;
    int background_mu_col;
    int background_mv_col;
    

};

void initialize(ggl_data_info *data_info)

#endif
