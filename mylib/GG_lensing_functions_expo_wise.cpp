#include<GG_lensing_functions_expo_wise.h>


void ggl_initialize(ggl_data_info *data_info)
{
    data_info->foreground_ra_col = 0;
    data_info->foreground_dec_col = 1;
    data_info->foreground_cos_dec_col = 2;
    data_info->foreground_z_col = 3;
    data_info->foreground_com_dist_col = 4;

    data_info->background_mg1_col = 0;
    data_info->background_mg2_col = 1;
    data_info->background_mn_col = 2;
    data_info->background_mu_col = 3;
    data_info->background_mv_col = 4;

    data_info->background_ra_col = 5;
    data_info->background_dec_col = 6;
    data_info->background_cos_dec_col = 7;
    data_info->background_z_col = 8;
    data_info->background_com_dist_col = 9;



    data_info->back_dz = 0.05;

    data_info->pos_inform_num = 6;

    data_info->gt_guess = new MY_FLOAT[data_info->pdf_guess_num];
    data_info->delta_sigma_guess = new MY_FLOAT[data_info->pdf_guess_num];
    data_info->mg_bin = new MY_FLOAT[data_info->mg_bin_num];

    data_info->sepration_bin = new MY_FLOAT[data_info->signal_pts_num+1];

    data_info->signal_chi_len = data_info->signal_pts_num*data_info->pdf_guess_num*data_info->mg_bin_num;
    data_info->chi_len = data_info->signal_chi_len*4 + data_info->signal_pts_num*3;
    // all pair data will be stacked into the last block in total_chi_count for the estimation of signals
    data_info->total_chi_count = new double[data_info->chi_len*(data_info->jack_num+1)];
    data_info->indi_chi_count = new double[data_info->chi_len];

    

}