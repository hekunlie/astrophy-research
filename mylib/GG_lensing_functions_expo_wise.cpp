#include<GG_lensing_functions_expo_wise.h>

void ggl_initialize(ggl_data_info *data_info)
{
    data_info->len_ra_col = 0;
    data_info->len_dec_col = 1;
    data_info->len_cos_dec_col = 2;
    data_info->len_z_col = 3;
    data_info->len_com_dist_col = 4;
    data_info->len_jackid_col = 5;
    
    data_info->src_mg1_col = 0;
    data_info->src_mg2_col = 1;
    data_info->src_mn_col = 2;
    data_info->src_mu_col = 3;
    data_info->src_mv_col = 4;

    data_info->src_ra_col = 5;
    data_info->src_dec_col = 6;
    data_info->src_cos_dec_col = 7;
    data_info->src_z_col = 8;


    // no len/src data array exists in memory
    data_info->len_expo_read_tag = 0;
    data_info->src_expo_read_tag = 0;


    data_info->back_dz = 0.05;

    data_info->pos_inform_num = 6;

    // when read a len exposure, it will calculate the separation between the 
    // len and src and labelled the src exposures.
    data_info->src_expo_needed_tag = new int[data_info->src_expo_num];
    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num,-1);

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


void ggl_read_len_exp(ggl_data_info *data_info, int len_expo_label)
{
    if(data_info->len_expo_read_tag == 1)
    {delete[] data_info->len_expo_data;}

    char set_name[50];
    int size;
    sprintf(set_name,"/data");
    
    size = data_info->len_data_col*data_info->len_data_row[len_expo_label];
    data_info->len_expo_data = new MY_FLOAT[size];
    read_h5(data_info->len_expo_path[len_expo_label], set_name, data_info->len_expo_data);
    
    data_info->len_expo_read_tag = 1;

#ifdef GGL_PROP_DIST_STACK
    // stack the signal in the proper distance coordinate
    int i, iz, ir;
    for(i=0; i<data_info->len_data_row[len_expo_label];i++)
    {   
        iz = i*data_info->len_data_col + data_info->len_z_col;
        ir = i*data_info->len_data_col + data_info->len_com_dist_col;
        data_info->len_expo_data[ir] = data_info->len_expo_data[ir]/(1+data_info->len_expo_data[iz]);
    }
#endif
}


void ggl_read_src_exp(ggl_data_info *data_info, int src_expo_label)
{
    if(data_info->src_expo_read_tag == 1)
    {delete[] data_info->src_expo_data;}

    char set_name[50];
    int size;
    sprintf(set_name,"/data");
    
    size = data_info->src_data_col*data_info->src_data_row[src_expo_label];
    data_info->src_expo_data = new MY_FLOAT[size];
    read_h5(data_info->src_expo_path[src_expo_label], set_name, data_info->src_expo_data);
    
    data_info->src_expo_read_tag = 1;
}

void find_src_needed(ggl_data_info *data_info, int len_expo_label)
{
    int fg, bkg;
    MY_FLOAT max_sep_theta;
    MY_FLOAT dra, ddec, sep_theta;
    MY_FLOAT len_ra, len_dec, len_cos_dec;

    initialize_arr(data_info->src_expo_needed_tag, data_info->src_expo_num, -1);

    for(fg=0; fg<data_info->len_data_row[len_expo_label]; fg++)
    {
        len_ra = data_info->len_expo_data[fg*data_info->len_data_col + data_info->len_ra_col];
        len_dec = data_info->len_expo_data[fg*data_info->len_data_col + data_info->len_dec_col];
        len_cos_dec = data_info->len_expo_data[fg*data_info->len_data_col + data_info->len_cos_dec_col];
        
        max_sep_theta = data_info->separation_bin[data_info->sep_bin_num+1]/
        data_info->len_expo_data[fg*data_info->len_data_col + data_info->len_com_dist_col]*1.2;

        for(bkg=0; bkg<data_info->src_expo_num; bkg++)
        {
            dra = (len_ra - data_info->src_pos_informs[bkg][0])*len_cos_dec;
            ddec = len_dec - data_info->src_pos_informs[bkg][1];
            sep_theta = (sqrt(dra*dra + ddec*ddec) - data_info->src_pos_informs[bkg][2])/180*Pi;

            if(sep_theta <= max_sep_theta){data_info->src_expo_needed_tag[bkg] = 1;}
        }
    }

}
void ggl_find_pair(ggl_data_info *data_info, int len_expo_label)
{
    int ibkg, bkg, ifg;
    MY_FLOAT len_ra, len_dec, src_ra, src_dec;
    MY_FLOAT sep_dist, sep_theta;
    int ir, bin_tag;
    
    for(bkg=0; bkg<data_info->src_expo_num; bkg++)
    {   
        if(data_info->src_expo_needed_tag[bkg]< 0){continue;}
        
        ggl_read_src_exp(data_info, bkg);

        for(ifg=0; ifg<data_info->len_data_row[len_expo_label]; ifg++)
        {   
            len_ra = data_info->len_expo_data[ifg*data_info->len_data_col + data_info->len_ra_col];
            len_dec = data_info->len_expo_data[ifg*data_info->len_data_col + data_info->len_dec_col];

            for(ibkg=0; ibkg<data_info->src_data_row[bkg]; ibkg++)
            {
                src_ra = data_info->src_expo_data[ibkg*data_info->src_data_col + data_info->src_ra_col];
                src_dec = data_info->src_expo_data[ibkg*data_info->src_data_col + data_info->src_dec_col];

                separation(len_ra, len_dec, src_ra, src_dec, sep_theta);
                sep_dist = sep_theta*data_info->len_expo_data[ifg*data_info->len_data_col + data_info->len_com_dist_col];

                bin_tag = -1;
                for(ir=0; ir<data_info->sep_bin_num; ir++)
                {
                    if(sep_dist >= data_info->separation_bin[ir] and sep_dist< data_info->separation_bin[ir+1])
                    {
                        bin_tag = ir;
                        break;
                    }
                }
                
                if(bin_tag > -1)
                {


                }
                

            }
        }
    }

}
