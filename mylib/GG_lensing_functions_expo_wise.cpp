#include<GG_lensing_functions_expo_wise.h>


void initialize(ggl_data_info *data_info)
{
    data_info->foreground_ra_col = 0;
    data_info->foreground_dec_col = 1;
    data_info->foreground_cos_dec_col = 2;
    data_info->foreground_z_col = 3;
    data_info->foreground_com_dist_col = 4;

    
    data_info->background_ra_col = 0;
    data_info->background_dec_col = 1;
    data_info->background_cos_dec_col = 2;
    data_info->background_z_col = 3;
    data_info->background_com_dist_col = 4;

    data_info->background_mg1_col = 5;
    data_info->background_mg2_col = 6;
    data_info->background_mn_col = 7;
    data_info->background_mu_col = 8;
    data_info->background_mv_col = 9;

    data_info->signal_chi_len = data_info->signal_pts_num;
    data_info->chi_len = data_info->signal_chi_len*4 + data_info->signal_pts_chi_len*3;

    data_info->total_chi_count = new double[data_info->chi_len*data_info->jack_num];
    data_info->indi_chi_count = new double[data_info->chi_len];

}